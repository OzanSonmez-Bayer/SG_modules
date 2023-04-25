# This script to calculate GCA and degressed EBV

library(plyr)
library(tidyverse)
library(asreml)
library(asremlPlus)

dat <- read.csv('Pepper/pheno.csv')

# Using wall thickness as an example for GBLUP calculation

########################################## Clean data ##############################

namelist <- paste('fruit',levels(interaction(c(1:9), c(1:3))), sep = '.')
namelist <- namelist[-length(namelist)]

dat_cln <- dat %>% 
  select(plot, Pedigree, Female, Male, namelist) %>% 
  mutate(avg_thickness = apply(dat[,namelist], 
                               MARGIN = 1, 
                               FUN = function(x){
                                 sum(x, na.rm = T)/sum(!(is.na(x) | x == 0))
                               })) %>% 
  select(plot, Pedigree, Female, Male, avg_thickness)


######################################### Calculate GCA ##############################

pd <- dat_cln %>% 
  select(Pedigree, Female, Male) %>% 
  unique()
pd <- pd[!duplicated(pd[,1]),]
names(pd) <- c('Pedigree', 'Female', 'Male')
pd$Selfing <- rep(0, rep = nrow(pd))

ped_matrix <- asreml.Ainverse(pd, fgen = c('Selfing',5))

mdl_blup <- asreml(fixed = avg_thickness ~ 1, 
                   random = ~ped(Pedigree), 
                   ginverse = list(Pedigree = ped_matrix$ginv),
                   na.method.Y = 'include',
                   na.method.X = 'include',
                   data = dat_cln)
sum_tab <- summary(mdl_blup)$varcomp

pred_ebv <- as.data.frame(predict(mdl_blup, classify = 'ped(Pedigree)', max.iter = 1)$predictions$pval)
pred_ebv <- pred_ebv %>% 
  mutate(mdl = NA, 
         mdl = if_else(Pedigree %in% unique(dat_cln$Pedigree), 'SCA', 'GCA'))

rnd_df <- data.frame(
                     BLUP = coef(mdl_blup)$random, 
                     BLUP.se = mdl_blup$vcoeff$random
                     )
colnames(rnd_df) <- c('BLUP', 'BLUP.se')
rnd_df$Pedigree <- gsub('ped(Pedigree)_', '', rownames(rnd_df), fixed = TRUE)

PEV <- rnd_df[,2]*sum_tab['R!variance',2]
reliability <- round(1-PEV/((1+ped_matrix$inbreeding)*sum_tab["ped(Pedigree)!ped",2]),4)

rnd_df$reliability <- reliability

pred_df <- pred_ebv %>% 
  join(rnd_df, by = 'Pedigree') %>% 
  mutate(dEBV = BLUP/reliability) %>% 
  filter(mdl == 'GCA')

save(pred_df, file = 'dEBV.Rdata')

############################################# GBLUP #############################
load('Pepper/dEBV.Rdata')
pheno_dat <- pred_df %>% 
  dplyr::select(Pedigree, dEBV) %>% 
  mutate(Pedigree = gsub('-', '.', Pedigree))
ped_list <- pheno_dat$Pedigree
pheno_dat <- data.frame(thk_dEBV = pheno_dat$dEBV)
rownames(pheno_dat) <- ped_list

# Step 3: Match Pedigrees from geno with pheno

identical(sort(rownames(geno_cln)), sort(rownames(pheno_dat)))

geno_sub <- geno_cln[rownames(geno_cln) %in% rownames(pheno_dat),]
geno_sub <- genElemnet(geno_sub) # replace markers with more than 2 alleles, i.e., TAATAA

map_sub <- map_dat %>% 
  select(X, CHROMOSOME, POSITION) %>% 
  filter(X %in% colnames(geno_sub)) %>% 
  mutate(POSITION = as.numeric(POSITION),
         CHROMOSOME = as.numeric(CHROMOSOME))
rownames(map_sub) <- map_sub$X
map_sub <- map_sub[,-1]
colnames(map_sub) <- c('chr', 'pos')


# Step 4: Create gpData object

gp <- create.gpData(geno = geno_sub, pheno = pheno_dat, pedigree = NULL, map = map_sub)
summary(gp)

# Step 5: Recode Loci and Impute Missing Genotypes
gp_num <- codeGeno(gp, label.heter = 'alleleCoding', maf = 0.05,
                   nmiss = 0.25, 
                   impute = TRUE, 
                   impute.type = 'random',
                   verbose = TRUE)
summary(gp_num)

# Step 6: Data Summary and Visualization
# library(ggplot2)
# library(RColorBrewer)
# library(LDheatmap)
# library(Hmisc)
# library(GGally)

# plotGenMap(gp_num, dense = F, nMarker = T)

# Step 6: Generate Kinship Matrix
addititve <- kin(gp_num, ret = 'add') # pedigree is required

realized <- kin(gp_num, ret = 'realized') # G matrix
summary(realized)
# dimension                     53 x 53 
# rank                          52 
# range of off-diagonal values  -0.5461042 -- 1.612514 
# mean off-diagonal values      -0.03773814 
# range of diagonal values      1.536 -- 3.135 
# mean diagonal values          1.962383

# Negative value at off-diagonal means that pairs of individuals that are less related than expected by random chance. 
# Rnak is 52, smaller than the number of genotypes, which indicates the existence of singularity. 
plot(realized, main = 'G matrix')

# Treat singularity in G matrix
library(Matrix)
realiszedPD <- nearPD(realized, keepDiag = T)
G <- matrix(realiszedPD[[1]]@x, nrow = realiszedPD[[1]]@Dim[1])
G <- G + diag(0.01, nrow(G))
attr(G, 'dimnames') <- realiszedPD[[1]]@Dimnames
class(G) <- 'relationshipMatrix'
str(G)
summary(G)
summary(eigen(G)$values)

plot(G, main = 'G Matrix')

# Save G matrix in the format that is comptable with Asreml
Gtab <- write.relationshipMatrix(G, sorting = 'ASReml', type = 'ginv')


# Step 7: Calculate GBLUP
# synbreed

Gblup <- gpMod(gp_num, model = 'BLUP', kin = G)
summary(Gblup)

cor(gp_num$pheno, Gblup$g, method = 'spearman') # 0.968
plot(Gblup$y, Gblup$g, pch = 20, cex = 1.6, col = 'gray', xlab = 'Observed', ylab = 'Predicted')
legend('topleft', paste('r_G = ', round(cor(gp_num$pheno, Gblup$g, method = 'spearman'), 3), sep = ''),
       col = 'red3',
       bty = 'n')

# Asreml
library(asreml)

pheno_dat_reml <- pheno_dat %>% 
  mutate(pedigree = rownames(pheno_dat))

gblup_reml <- asreml(fixed = thk_dEBV ~ 1, random = ~pedigree,
                     ginverse=list(pedigree = Gtab),
                     data = pheno_dat_reml)

reml_blup <- as.data.frame(predict(gblup_reml, classify = 'pedigree', maxiter = 1)$predictions$pval)

cor(reml_blup$predicted.value, Gblup$g, method = 'spearman')
plot(gp_num$pheno, Gblup$g, pch = 20, cex = 1.6, col = 'gray', xlab = 'Observed', ylab = 'Predicted', main = 'Prediction from Synbreed')
legend('topleft', paste('r_synbreed = ', round(cor(gp_num$pheno, Gblup$g, method = 'spearman'), 3), sep = ''),
       col = 'red3',
       bty = 'n')

plot(gp_num$pheno, reml_blup$predicted.value, col = 'blue', pch = 20, cex = 1.6,  xlab = 'Observed', ylab = 'Predicted', 
     main = 'Prediction from Asreml')
legend('topleft', paste('R_asreml = ', round(cor(gp_num$pheno, reml_blup$predicted.value, method = 'spearman'), 3), sep = ''),
       col = 'red',
       bty = 'n')



################################# Test GBLUP with New Lines ##################################
library(synbreed)
library(readxl)
source('genoFunctions.R')
new_pheno <- read_excel('Pepper/20190710 Fruit wall thickness measurements.xlsx')
new_geno <- read.csv('Pepper/Pepper_geno.csv')

# Step 1: Clean pheno

pheno_lsm <- new_pheno %>% 
  select(`Ave. Wall Thick`, pedigree) %>% 
  group_by(pedigree) %>% 
  summarise(LSM_thk = mean(`Ave. Wall Thick`, na.rm = T)) %>% 
  filter(!is.na(pedigree)) %>% 
  mutate(pedigree = gsub('-', '.', pedigree))

# Step 2: Clean geno
chr_pos <- new_geno %>% 
  select(MarkerName, Chr, Pos)

geno_t <- t(new_geno[, 4:ncol(new_geno)])
colnames(geno_t) <- new_geno$MarkerName
geno_t <- genElemnet(geno_t) # remove elements with "|" and more than one alleles
geno_t[which(geno_t == '??', arr.ind = TRUE)] <- NA

identical(sort(as.character(pheno_lsm$pedigree)), sort(rownames(geno_t)))

# Step 3: Match Pedigrees from Geno with Pheno

geno_sub <- geno_t[rownames(geno_t) %in% pheno_lsm$pedigree,]
map_sub <- chr_pos %>% 
  select(MarkerName, Chr, Pos) %>% 
  filter(MarkerName %in% colnames(geno_sub) & !is.na(Chr)) %>% 
  mutate(Chr = as.numeric(Chr),
         Pos = as.numeric(Pos))
rownames(map_sub) <- map_sub$MarkerName
map_sub <- map_sub[,-1]
colnames(map_sub) <- c('chr', 'pos')

pheno_sub <- pheno_lsm %>% 
  filter(as.character(pedigree) %in% rownames(geno_sub)) %>% 
  as.data.frame()

ped_list <- pheno_sub$pedigree
pheno_sub <- data.frame(LSM_thk = pheno_sub$LSM_thk) # pedigree has to be row names. columns in pheno are  traits
rownames(pheno_sub) <- ped_list

identical(sort(rownames(pheno_sub)), sort(rownames(geno_sub)))

# Step 4: Create a Pedigree 
ped_dat <- data.frame(ID = rownames(pheno_sub)[!is.na(rownames(pheno_sub))],
                      Par1 = 0,
                      Par2 = 0)

# ped_dat[1,2:3] <- strsplit(as.character(ped_dat[1,1]), '\\+')[[1]] # hard coded line
ped_dat$gener <- rep(0, nrow(ped_dat))

ped_file <- create.pedigree(ped_dat$ID, ped_dat$Par1, ped_dat$Par2, ped_dat$gener)

plot(ped_file)


# Step 5: Create gpData object
gp <- create.gpData(geno = geno_sub, pheno = pheno_sub, pedigree = ped_file, map = NULL)
summary(gp)

# check the missing values for each marker
nmissing_ind <- apply(geno_sub, 2, function(x){sum(x == 'NA')/length(x)})
hist(nmissing_ind, xlab = 'Percentage of Missing Values', main = '')

# Step 6: Recode Loci and Impute Missing Genotypes
gp_num <- codeGeno(gp, label.heter = 'alleleCoding',
                   maf = 0.05,
                   nmiss = 0.25,
                   impute = TRUE, 
                   impute.type = 'random',
                   verbose = TRUE)

# Step 7: Generate A Matrix and G Matrix
Additive <- kin(gp_num, ret = 'add')
summary(Additive)

# dimension                     82 x 82 
# rank                          82 
# range of off-diagonal values  0 -- 0 
# mean off-diagonal values      0 
# range of diagonal values      1 -- 1 
# mean diagonal values          1

Realized <- kin(gp_num, ret = 'realized')
Realized <- kin(gp_num, ret = 'realizedAB') # weighted G
summary(Realized)

# dimension                     82 x 82 
# rank                          81 
# range of off-diagonal values  -0.6220297 -- 2.189381 
# mean off-diagonal values      -0.02350584 
# range of diagonal values      0.9219 -- 2.383 
# mean diagonal values          1.903973 

# par(mfrow = c(1,2))
plot(Realized, main = 'G Matrix')
plot(Additive, main = 'A Matrix')

# G matrix is not full-rank, which indicates that there is singularity. 
library(Matrix)
RealizedPD <- nearPD(Realized, keepDiag = T)
G <- matrix(RealizedPD[[1]]@x, nrow = RealizedPD[[1]]@Dim[1])
G <- G + diag(0.01, nrow(G))
attr(G, 'dimnames') <- RealizedPD[[1]]@Dimnames
class(G) <- 'relationshipMatrix'
str(G)
summary(G)

# dimension                     82 x 82 
# rank                          82 
# range of off-diagonal values  -0.6220297 -- 2.189381 
# mean off-diagonal values      -0.02350584 
# range of diagonal values      0.9319 -- 2.393 
# mean diagonal values          1.913973 

plot(G, main = 'G Matrix')

par(mfrow = c(1,1))
boxplot(as.vector(G)~as.vector(Additive), xlab = 'A', ylab = 'G')

# Save G matrix in the format that is comptable with Asreml
Gtab <- write.relationshipMatrix(G, sorting = 'ASReml', type = 'ginv')
Atab <- write.relationshipMatrix(Additive, sorting = 'ASReml', type = 'ginv')


# Step 8: A-BLUP and G-BLUP
Ablup <- gpMod(gp_num, model = 'BLUP', kin = Additive, trait = 'LSM_thk')
summary(Ablup)

# Linear Coefficients:
#              Estimate Std. Error
# (Intercept)    7.521      0.065
# 
# Variance Coefficients:
#       Estimate Std. Error
# kinTS    0.173      0.027
# In       0.173      0.027

Gblup <- gpMod(gp_num, model = 'BLUP', kin = G, trait = 1)
summary(Gblup)

# Linear Coefficients:
#              Estimate Std. Error
# (Intercept)    7.521       0.05
# 
# Variance Coefficients:
#         Estimate Std. Error
# kinTS    0.087      0.042
# In       0.205      0.056

# Step 9: Correlation for Accuracy
# Note g is raw blup
# predicted value = g + intercept from linear coefficients
# Compare Rank
cor(Ablup$y, Ablup$g,method = 'spearman') # 1
cor(Gblup$y, Gblup$g, method = 'spearman') # 0.857

# Compare actual predicted value vs LSM
cor(Ablup$y, Ablup$g+7.521) # 1
cor(Gblup$y, Gblup$g+7.521) # 0.882

plot(Gblup$g+7.521, Gblup$y, pch = 20, cex = 1.6, col = 'gray',
     xlab = 'Observed', 
     ylab = 'Predicted')
points(Ablup$g+7.521, Ablup$y, pch = 1, cex = 0.6, col = 'blue',
       xlab = 'Observed', 
       ylab = 'Predicted')
legend('topleft', paste('r_G = ', round(cor(Gblup$y, Gblup$g+7.521), 3), sep = ''),
       col = 'red3',
       bty = 'n')
legend('bottomright', paste('r_A = ', round(cor(Ablup$y, Ablup$g+7.521), 3), sep = ''),
       col = 'red',
       bty = 'n')

# ABLUP has higher correlation that GBLUP with LSM, it indicates that ABLUP is more similar to the original observations. 
# Correlation may not be adequate evaluation of the predictive ability of the models


############################ Cross-Validation #################

# CV GBLUP
# Create a K-Fold Cross-Validation
k_fold <- 5

set.seed(20191023)
new_ph <- as.data.frame(gp_num$pheno) # start by from gpData object after QC
# new_ph$pedigree <- rownames(gp_num$pheno)
new_ph$fold <- sample(1:k_fold, size = nrow(new_ph), replace = T)
new_y <- as.matrix(new_ph[, 1, drop = F]) %*% t(rep(1,k_fold))
colnames(new_y) = paste0('Y', 1:k_fold)
new_ph <- cbind(new_ph, new_y)

missingInFold <- function(r){
  fold = r['fold']
  index = fold + 2
  r[index] = NA
  return(r)
}

new_ph_2 <- t(apply(new_ph, 1, missingInFold))
new_ph_df <- data.frame(new_ph_2)
new_ph_df$pedigree <- factor(rownames(new_ph_df))

G.giv <- write.relationshipMatrix(G, file = NULL,
                                  sorting = c('ASReml'), 
                                  type = c('ginv'), digits = 10)

head(attr(G.giv, 'rowNames'))
names(G.giv) <- c('row', 'column', 'coefficient')
head(G.giv)

# The following set up could be used for multiple traits
result.list <- list()
pred.list <- list()
for (trait in names(new_ph_df) [3:(3+k_fold-1)]){
  asr <- asreml(fixed = get(trait) ~ 1,
                random = ~giv(pedigree, var = T),
                ginverse = list(pedigree = G.giv),
                rcov = ~units,
                data = new_ph_df,
                na.method.X = 'include')
  result.list[[trait]] <- asr
  pred.list[[trait]] <- summary(asr, all = T)$coef.rand
  rm(list = 'asr')
}

# Summarize cross validation
get_varcomps <- function(c){
  vcs <- summary(c)$varcomp$component
  names(vcs) <- c('Vg', 'Vresid')
  return (vcs)
}

varcomps <- sapply(result.list, FUN = get_varcomps)
varcomps <- as.data.frame(t(varcomps))
varcomps$fold <- as.numeric(gsub('Y', '', rownames(varcomps)))

# get the mean for each estimation set for each fold to use as mean of predictions
mu <- colMeans(new_ph_df[, paste('Y', 1:k_fold, sep = "")], na.rm = T)
mu <- as.data.frame(mu)
mu$fold <- as.numeric(gsub('Y', '', rownames(mu)))

vc_mu <- merge(varcomps, mu)

# for each fold, get the predictions only for the held out observations
get_pred <- function(r){
  ID <- as.character(r['pedigree'])
  ID2 <- paste0('giv(pedigree, var = T)_', ID)
  trait <- paste0("Y", r['fold'])
  return(pred.list[[trait]][ID2, ])
}

preds <- t(apply(new_ph_df, 1, get_pred))

pred_obs <- merge(new_ph_df[,c('LSM_thk.1', 'fold')],
                   as.data.frame(preds), by = 0)

pred_obs <- merge(pred_obs, vc_mu, by = 'fold')
pred_obs$pred <- pred_obs$mu + pred_obs$solution
colnames(pred_obs) <- c('fold', 'Pedigree', 'thk_observed', 'BLUP', 'BLUP.se', 'z', 'Vg', 'Vresid', 'mu', 'thk_pred')
pred_obs$reliability <- with(pred_obs, 1 - (BLUP.se^2*Vresid)/(Vg * (1+mean(diag(G))))) 
pred_obs$h2 <- with(pred_obs, Vg/(Vg+Vresid))

# Alternative way to get predicted value
# predict(result.list[[1]], classify = 'giv(pedigree, var = T)')$predictions$pval

# Summary Statistics
# PredCor <- with(pred_obs, cor(thk_observed, thk_pred))
# RankCor <- with(pred_obs, cor(thk_observed, thk_pred, method = 'spearman'))
# Bias <- coef(lm(pred_obs$thk_observed ~ pred_obs$thk_pred))[2]
# MSE <- mean((pred_obs$thk_observed - pred_obs$thk_pred)**2)
# pred_sum <- rbind(PredCor, RankCor, Bias, MSE)
# pred_sum

# Comments: double check correlation

pred_sum <- pred_obs %>%
  group_by(fold) %>%
  mutate(PredCor = cor(thk_observed, thk_pred),
         RankCor = cor(thk_observed, thk_pred, method = 'spearman'),
         MSE = mean(thk_observed - thk_pred)**2,
         Bias = coef(lm(thk_observed ~ thk_pred))[2]) %>% 
  ungroup() %>% 
  select(PredCor, RankCor, MSE, Bias)
colMeans(pred_sum)
#   PredCor    RankCor        MSE       Bias 
# 0.29599149 0.20618222 0.03882574 1.69511635  

# Check the distribution of observed and predicted value
thk_obs <- data.frame(avg_thk = pred_obs$thk_observed)
thk_pred <- data.frame(avg_thk = pred_obs$thk_pred)
thk_obs$type <- 'Observed'
thk_pred$type <- 'GBLUP'
thk_df <- rbind(thk_obs, thk_pred)
ggplot(thk_df, aes(avg_thk, fill = type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity',bins = 10)

ggplot(thk_df, aes(avg_thk, fill = type)) + 
  geom_density(alpha = 0.2)

# Bias is 1.69 indicates that 'deflation' meaning less variance among the predicted than observed values. 
# Overlapped density plot is congurent with the result 


############################## GBLUP Combining Old and New Population ###################################
load('Pepper/dEBV.Rdata')
geno_new <- read.csv('Pepper/Pepper_geno.csv')
geno_old <- read.csv('Pepper/Pepper_geno_old.csv')

# Clean Genotype Files
# remove duplicated genotypes
geno_old <- geno_old %>% 
  select(colnames(geno_old)[!colnames(geno_old) %in% intersect(colnames(geno_new), colnames(geno_old))[-1]]) 

# Combine old and new populations
geno_file <- merge(geno_new, geno_old, by = 'MarkerName')

chr_pos <- geno_file %>% 
  select(MarkerName, Chr, Pos)

geno_t <- t(geno_file[, 4:ncol(geno_file)])
colnames(geno_t) <- geno_file$MarkerName
geno_t <- genElemnet(geno_t)
geno_t[which(geno_t == '??', arr.ind = TRUE)] <- NA


# Clean Pheno Files

colnames(pheno_lsm) <- c('Pedigree', 'Avg_thk')
pheno_old <- pred_df %>% 
  select(Pedigree, predicted.value) %>% 
  mutate(Pedigree = gsub('-', '.', Pedigree))

colnames(pheno_old) <- colnames(pheno_lsm)
pheno_file <- as.data.frame(rbind(pheno_lsm, pheno_old))
pheno_file <- pheno_file[!duplicated(pheno_file$Pedigree),]
row_nm <- pheno_file$Pedigree
pheno_file <- data.frame(Avg_thk = pheno_file$Avg_thk)
rownames(pheno_file) <- row_nm

# Match Pedigrees from Geno with Pheno
geno_sub <- geno_t[rownames(geno_t) %in% rownames(pheno_file),]
map_sub <- chr_pos %>% 
  select(MarkerName, Chr, Pos) %>% 
  filter(MarkerName %in% colnames(geno_sub) & !is.na(Chr)) %>% 
  mutate(Chr = as.numeric(Chr),
         Pos = as.numeric(Pos))
rownames(map_sub) <- map_sub$MarkerName
map_sub <- map_sub[,-1]
colnames(map_sub) <- c('chr', 'pos')
pheno_sub <- pheno_file %>% 
  mutate(pedigree = rownames(pheno_file)) %>% 
  filter(pedigree %in% rownames(geno_sub))

ped_list <- pheno_sub$pedigree
pheno_sub <- pheno_sub %>% 
  select(Avg_thk) %>% 
  as.data.frame()
rownames(pheno_sub) <- ped_list

identical(sort(rownames(pheno_sub)), sort(rownames(geno_sub)))

# Create gpData object
gp <- create.gpData(geno = geno_sub, pheno = pheno_sub, pedigree = NULL, map = map_sub)
summary(gp)
plot(gp$map)

# Impute Missing Data
gp_num <- codeGeno(gp, label.heter = 'alleleCoding',
                   maf = 0.05,
                   nmiss = 0.25,
                   impute = TRUE,
                   impute.type = 'random',
                   verbose = TRUE)
# step 1  : 521 marker(s) removed with > 25 % missing values 
# step 2  : Recoding alleles 
# step 4  : 2614 marker(s) removed with maf < 0.05 
# step 7  : Imputing of missing values 
# step 7d : Random imputing of missing values 
# approximate run time for random imputation  0.297  seconds 
# step 8  : Recode alleles due to imputation 
# step 9  : 1 marker(s) removed with maf < 0.05 
# step 10 : No duplicated markers removed 
# End     : 1074 marker(s) remain after the check
# 
# Summary of imputation 
# total number of missing values                : 6319 
# number of random imputations                  : 6319 


# Generate G Matrix
realized <- kin(gp_num, ret = 'realized')
summary(realized)

# Remove singularity in the matrix
library(Matrix)

realizedPD <- nearPD(realized, keepDiag = T)
G <- matrix(realizedPD[[1]]@x, nrow = realizedPD[[1]]@Dim[1])
G <- G + diag(0.01, nrow(G))
attr(G, 'dimnames') <- realizedPD[[1]]@Dimnames
class(G) <- 'relationshipMatrix'
str(G)
summary(G)
# dimension                     125 x 125 
# rank                          125 
# range of off-diagonal values  -0.5519344 -- 2.024316 
# mean off-diagonal values      -0.01548031 
# range of diagonal values      1.01 -- 2.804 
# mean diagonal values          1.929558 
plot(G)

Gtab <- write.relationshipMatrix(G, sorting = 'ASReml', type = 'ginv')

# Cross Validation
k_fold <- 5

set.seed(20191028)
new_ph <- as.data.frame(gp_num$pheno) # start by from gpData object after QC
# new_ph$pedigree <- rownames(gp_num$pheno)
new_ph$fold <- sample(1:k_fold, size = nrow(new_ph), replace = T)
new_y <- as.matrix(new_ph[, 1, drop = F]) %*% t(rep(1,k_fold))
colnames(new_y) = paste0('Y', 1:k_fold)
new_ph <- cbind(new_ph, new_y)

new_ph_2 <- t(apply(new_ph, 1, missingInFold))
new_ph_df <- data.frame(new_ph_2)
new_ph_df$pedigree <- factor(rownames(new_ph_df))

G.giv <- write.relationshipMatrix(G, file = NULL,
                                  sorting = c('ASReml'), 
                                  type = c('ginv'), digits = 10)

head(attr(G.giv, 'rowNames'))
names(G.giv) <- c('row', 'column', 'coefficient')
head(G.giv)

result.list <- list()
pred.list <- list()
for (trait in names(new_ph_df) [3:(3+k_fold-1)]){
  asr <- asreml(fixed = get(trait) ~ 1,
                random = ~giv(pedigree, var = T),
                ginverse = list(pedigree = G.giv),
                rcov = ~units,
                data = new_ph_df,
                na.method.X = 'include')
  result.list[[trait]] <- asr
  pred.list[[trait]] <- summary(asr, all = T)$coef.rand
  rm(list = 'asr')
}

varcomps <- sapply(result.list, FUN = get_varcomps)
varcomps <- as.data.frame(t(varcomps))
varcomps$fold <- as.numeric(gsub('Y', '', rownames(varcomps)))

# get the mean for each estimation set for each fold to use as mean of predictions
mu <- colMeans(new_ph_df[, paste('Y', 1:k_fold, sep = "")], na.rm = T)
mu <- as.data.frame(mu)
mu$fold <- as.numeric(gsub('Y', '', rownames(mu)))

vc_mu <- merge(varcomps, mu)

preds <- t(apply(new_ph_df, 1, get_pred))

pred_obs <- merge(new_ph_df[,c('Avg_thk.1', 'fold')],
                  as.data.frame(preds), by = 0)

pred_obs <- merge(pred_obs, vc_mu, by = 'fold')
pred_obs$pred <- pred_obs$mu + pred_obs$solution
colnames(pred_obs) <- c('fold', 'Pedigree', 'thk_observed', 'BLUP', 'BLUP.se', 'z', 'Vg', 'Vresid', 'mu', 'thk_pred')
pred_obs$reliability <- with(pred_obs, 1 - (BLUP.se^2*Vresid)/(Vg * (1+mean(diag(G))))) 
pred_obs$h2 <- with(pred_obs, Vg/(Vg+Vresid))

pred_sum <- pred_obs %>%
  group_by(fold) %>%
  mutate(PredCor = cor(thk_observed, thk_pred),
         RankCor = cor(thk_observed, thk_pred, method = 'spearman'),
         MSE = mean(thk_observed - thk_pred)**2,
         Bias = coef(lm(thk_observed ~ thk_pred))[2]) %>% 
  ungroup() %>% 
  select(PredCor, RankCor, MSE, Bias)
colMeans(pred_sum)
#    PredCor     RankCor         MSE        Bias 
# 0.222009361 0.173139261 0.009539906 0.764058268 

thk_obs <- data.frame(avg_thk = pred_obs$thk_observed)
thk_pred <- data.frame(avg_thk = pred_obs$thk_pred)
thk_obs$type <- 'Observed'
thk_pred$type <- 'GBLUP'
thk_df <- rbind(thk_obs, thk_pred)
ggplot(thk_df, aes(avg_thk, fill = type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity',bins = 10)

ggplot(thk_df, aes(avg_thk, fill = type)) + 
  geom_density(alpha = 0.2)

############################### 90/10 Cross Validation #########################

set.seed(20191028)
training <- sample(1:nrow(gp_num$pheno), round(nrow(gp_num$pheno)*0.9,0), replace = F)
new_ph <- as.data.frame(gp_num$pheno) 
new_ph$thk_cv <- new_ph$Avg_thk.1
new_ph$thk_cv[-training] <- NA
new_ph$pedigree <- rownames(new_ph)

G.giv <- write.relationshipMatrix(G, file = NULL,
                                  sorting = c('ASReml'), 
                                  type = c('ginv'), digits = 10)

head(attr(G.giv, 'rowNames'))
names(G.giv) <- c('row', 'column', 'coefficient')
head(G.giv)

mld_gblup <- asreml(fixed = thk_cv ~ 1,
                    random = ~giv(pedigree, var = T),
                    ginverse = list(pedigree = G.giv),
                    rcov = ~units,
                    data = new_ph,
                    na.method.X = 'include')
summary(mld_gblup)$varcomp # 0.101/(0.101 + 0.171) = 0.37

pred_gebv <- predict(mld_gblup, classify = 'giv(pedigree, var = T)')$predictions$pval

pred_observed <- merge(pred_gebv, new_ph, by = 'pedigree')

# All the obs
cor(pred_observed$predicted.value, pred_observed$Avg_thk.1) # 0.84

# Validation 
cor(pred_observed$predicted.value[is.na(pred_observed$thk_cv)], 
    pred_observed$Avg_thk.1[is.na(pred_observed$thk_cv)], 
    method = 'spearman') # 0.18
