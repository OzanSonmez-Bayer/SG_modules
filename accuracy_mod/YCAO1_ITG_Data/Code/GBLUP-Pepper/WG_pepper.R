# Cross validation with old and new lines
# 02/25/2020

library(snpReady)
library(tidyverse)
library(BGLR)
library(MASS)
library(synbreed)
library(readxl)
library(Matrix)
library(ggcorrplot)

source('Code/Shared/genoFunctions.R')
##################### Load Phenotypic Data #################

load('Pepper/dEBV.Rdata')
pheno_old <- pred_df
pheno_new <- read_excel('Pepper/20190710 Fruit wall thickness measurements.xlsx')

# LSM from new lines
pheno_new <- pheno_new %>% 
  mutate(pedigree = gsub('-', '\\.', pedigree)) %>% 
  group_by(pedigree) %>% 
  summarise(thickness_lsm = mean(`Ave. Wall Thick`, na.rm = T))

# Combine old and new lines

pheno <- pheno_old %>% 
  dplyr::select(Pedigree, predicted.value) %>% 
  rename(pedigree  = Pedigree, 
         thickness_lsm = predicted.value) %>% 
  dplyr::bind_rows(pheno_new) %>% 
  mutate(pedigree = gsub('-', '\\.', pedigree))

pheno <- pheno[!duplicated(pheno$pedigree), ]
pheno <- pheno[which(!is.na(pheno$pedigree)), ]
pheno_df <- data.frame(LSM = pheno$thickness_lsm)
rownames(pheno_df) <- pheno$pedigree


################### Load Genomic Data ####################

geno_old <- read.csv('Pepper/Pepper_geno_old.csv')
geno_new <- read.csv('Pepper/Pepper_geno.csv')  

geno_old <- geno_old %>% 
  dplyr::select(colnames(geno_old)[!colnames(geno_old) %in% intersect(colnames(geno_new), colnames(geno_old))[-1]]) 

# Combine old and new populations
geno_file <- merge(geno_new, geno_old, by = 'MarkerName')
geno_t <- t(geno_file[, 4:ncol(geno_file)])
colnames(geno_t) <- geno_file$MarkerName
geno_t <- genElemnet(geno_t)
geno_t[which(geno_t == '??', arr.ind = TRUE)] <- NA

################# Load Marker Map ########################

chr_pos <- geno_file %>% 
  dplyr::select(MarkerName, Chr, Pos)


################ Match geno with pheno ###################

geno_sub <- geno_t[rownames(geno_t) %in% rownames(pheno_df),]
map_sub <- chr_pos %>% 
  dplyr::select(MarkerName, Chr, Pos) %>% 
  filter(MarkerName %in% colnames(geno_sub) & !is.na(Chr)) %>% 
  mutate(Chr = as.numeric(Chr),
         Pos = as.numeric(Pos))
rownames(map_sub) <- map_sub$MarkerName
map_sub <- map_sub[,-1]
colnames(map_sub) <- c('chr', 'pos')
pheno_sub <- pheno_df %>% 
  mutate(pedigree = rownames(pheno_df)) %>% 
  filter(pedigree %in% rownames(geno_sub))

ped_list <- pheno_sub$pedigree
pheno_sub <- pheno_sub %>% 
  dplyr::select(LSM) %>% 
  as.data.frame()
rownames(pheno_sub) <- ped_list

identical(sort(rownames(pheno_sub)), sort(rownames(geno_sub)))

geno_sub[which(geno_sub == 'NA', arr.ind = TRUE)] <- NA

############## Generate Kinship Matrix ####################

##### Create gpData object
gp <- create.gpData(geno = geno_sub, pheno = pheno_sub, pedigree = NULL, map = map_sub)
summary(gp)
plot(gp$map)

##### Impute Missing Data
gp_num <- codeGeno(gp, label.heter = 'alleleCoding',
                   maf = 0.05,
                   nmiss = 0.25,
                   impute = TRUE,
                   impute.type = 'random',
                   verbose = TRUE)

##### Method 1: VanRaden, G_Base

realized <- kin(gp_num, ret = 'realized')
summary(realized)

# Remove singularity in the matrix

realizedPD <- nearPD(realized, keepDiag = T)
G_base <- matrix(realizedPD[[1]]@x, nrow = realizedPD[[1]]@Dim[1])
G_base <- G_base + diag(0.01, nrow(G_base))
attr(G_base, 'dimnames') <- realizedPD[[1]]@Dimnames
class(G_base) <- 'relationshipMatrix'
str(G_base)
summary(G_base)
# dimension                     125 x 125 
# rank                          125 
# range of off-diagonal values  -0.5662422 -- 2.028777 
# mean off-diagonal values      -0.01548537 
# range of diagonal values      0.9796846 -- 2.800697 
# mean diagonal values          1.930186 

plot(G_base)


##### Method 2: Weight is expected variance of each SNP, w = 2p_i(1-p_i)
G_snp <- weightedG_1(gp_num$geno)
# ped_name <- colnames(G_snp)
if (is.singular.matrix(G_snp)){
  G_snp_1 <- nearPD(G_snp, keepDiag = T)
  G_snp <- matrix(G_snp_1[[1]]@x, nrow = G_snp_1[[1]]@Dim[1])
  G_snp <- G_snp + diag(0.01, nrow(G_snp))
  
  attr(G_snp, 'dimnames') <- list(ped_name, ped_name)
  class(G_snp) <- 'relationshipMatrix'
}
# plot(G_snp)
ggcorrplot(G_snp[, ncol(G_snp):1]) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  ggtitle('Weighted Kinship Matrix G_SNP: Eexpected Variance')


##### Method 3: Bayesian LASSO 
G_BL <- weightedG_BL(phenodata = pheno_sub$LSM, markerM = gp_num$geno, 
                    r = 0.64, nIter = 20000, burnIn = 2000, outPath = 'Pepper/BL_')
ggcorrplot(G_BL[, ncol(G_BL):1]) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  ggtitle('Weighted Kinship Matrix G_BL: Bayesian LASSO')  


############## Cross Validation ####################
k_fold <- 5

set.seed(20200225)
new_ph <- as.data.frame(gp_num$pheno)
new_ph$fold <- sample(1:k_fold, size = nrow(new_ph), replace = T)
new_y <- as.matrix(new_ph[, 1, drop = F]) %*% t(rep(1, k_fold))
colnames(new_y) = paste0('Y', 1:k_fold)
new_ph <- cbind(new_ph, new_y)

new_ph_2 <- t(apply(new_ph, 1, missingInFold))
new_ph_df <- data.frame(new_ph_2)
new_ph_df$pedigree <- factor(rownames(new_ph_df))

G <- G_base # Switch to a different kinship matrix
G.giv_base <- write.relationshipMatrix(G, 
                                       file = NULL, 
                                       sorting = c('ASReml'),
                                       type = c('ginv'), 
                                       digits = 10)
head(attr(G.giv_base, 'rowNames'))
names(G.giv_base) <- c('row', 'column', 'coefficient')

result.list_base <- list()
pred.list_base <- list()
for (trait in names(new_ph_df)[3:(3 + k_fold - 1)]){
  asr <- asreml(fixed = get(trait) ~ 1,
                random = ~giv(pedigree, var = T), 
                ginverse = list(pedigree = G.giv_base), 
                rcov = ~units,
                data = new_ph_df, 
                na.method.X = 'include')
  result.list_base[[trait]] <- asr
  pred.list_base[[trait]] <- summary(asr, all = T)$coef.rand
  rm(list = 'asr')
}

varcomps <- sapply(result.list_base, FUN = get_varcomps)
varcomps <- as.data.frame(t(varcomps))
varcomps$fold <- as.numeric(gsub('Y', '', rownames(varcomps)))

mu <- colMeans(new_ph_df[, paste('Y',  1:k_fold, sep = '')], na.rm = T)
mu <- as.data.frame(mu)
mu$fold <- as.numeric(gsub('Y', '', rownames(mu)))

vc_mu <- merge(varcomps, mu)
preds <- t(apply(new_ph_df, 1, get_pred, pred.list_base))

pred_obs <- merge(new_ph_df[,c('LSM.1', 'fold')],
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
  dplyr::select(PredCor, RankCor, MSE, Bias)
colMeans(pred_sum)
#   PredCor     RankCor         MSE        Bias 
# 0.246744795 0.233867041 0.006860577 0.615554972 

thk_obs <- data.frame(avg_thk = pred_obs$thk_observed)
thk_pred <- data.frame(avg_thk = pred_obs$thk_pred)
thk_obs$type <- 'Observed'
thk_pred$type <- 'GBLUP'
thk_df <- rbind(thk_obs, thk_pred)
ggplot(thk_df, aes(avg_thk, fill = type)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity',bins = 10)

ggplot(thk_df, aes(avg_thk, fill = type)) + 
  geom_density(alpha = 0.2)



######### Different Marker Density ##################

ncol(geno_sub) # 4210 markers

# randomly sample certain percentage of markers from each Chromosome

set.seed(20200306)
k_fold <- 5

# 80% markers from each chromosome

map_80 <- sampleMark(dat = map_sub, perct = 0.8)


# Geno file
geno_80 <- geno_sub[, which(colnames(geno_sub) %in% rownames(map_80))]

G_80 <- unweightedGFun(geno_file = geno_80, pheno_file = pheno_sub, map_file = map_80, plot_draw = FALSE)

relatedness_corr(mat1 = G_80$G, mat2 = G_base)

# K-fold Cross Validation
pheno_80 <- kfoldValidation_pheno(k_fold = 5, pheno_dat = G_80$gp_obj$pheno)

G <- G_80$G # Switch to a different kinship matrix
G.giv_base <- write.relationshipMatrix(G, 
                                       file = NULL, 
                                       sorting = c('ASReml'),
                                       type = c('ginv'), 
                                       digits = 10)
head(attr(G.giv_base, 'rowNames'))
names(G.giv_base) <- c('row', 'column', 'coefficient')

result.list_base <- list()
pred.list_base <- list()
for (trait in names(pheno_80)[3:(3 + k_fold - 1)]){
  asr <- asreml(fixed = get(trait) ~ 1,
                random = ~giv(pedigree, var = T), 
                ginverse = list(pedigree = G.giv_base), 
                rcov = ~units,
                data = pheno_80, 
                na.method.X = 'include')
  result.list_base[[trait]] <- asr
  pred.list_base[[trait]] <- summary(asr, all = T)$coef.rand
  rm(list = 'asr')
}

pred_80 <- gblup_meta(pheno_80, G, result.list_base, pred.list_base)

accuracy_80 <- validationAccuracy(pred_80)
colMeans(accuracy_80)


# Test marker density from 10% to 100%

marker_density <- seq(0.1, 1, 0.1)
accuracy_df <- data.frame()
for (i in 1:length(marker_density)){
  print(i)
  
  set.seed(20200306)
  
  # 80% markers from each chromosome
  
  map_temp <- sampleMark(dat = map_sub, perct = marker_density[i])
  
  
  # Geno file
  geno_temp <- geno_sub[, which(colnames(geno_sub) %in% rownames(map_temp))]
  
  G_temp <- unweightedGFun(geno_file = geno_temp, pheno_file = pheno_sub, map_file = map_temp, plot_draw = FALSE)
  
  r <- relatedness_corr(mat1 = G_temp$G, mat2 = G_base)
  
  # K-fold Cross Validation
  pheno_temp <- kfoldValidation_pheno(k_fold = 5, pheno_dat = G_temp$gp_obj$pheno)
  
  G <- G_temp$G # Switch to a different kinship matrix
  G.giv_base <- write.relationshipMatrix(G, 
                                         file = NULL, 
                                         sorting = c('ASReml'),
                                         type = c('ginv'), 
                                         digits = 10)
  head(attr(G.giv_base, 'rowNames'))
  names(G.giv_base) <- c('row', 'column', 'coefficient')
  
  result.list_base <- list()
  pred.list_base <- list()
  for (trait in names(pheno_temp)[3:(3 + k_fold - 1)]){
    asr <- asreml(fixed = get(trait) ~ 1,
                  random = ~giv(pedigree, var = T), 
                  ginverse = list(pedigree = G.giv_base), 
                  rcov = ~units,
                  data = pheno_temp, 
                  na.method.X = 'include')
    result.list_base[[trait]] <- asr
    pred.list_base[[trait]] <- summary(asr, all = T)$coef.rand
    rm(list = 'asr')
  }
  
  pred_temp <- gblup_meta(pheno_temp, G, result.list_base, pred.list_base)
  
  accuracy_temp <- validationAccuracy(pred_temp)
  acc <- data.frame(t(colMeans(accuracy_temp)))
  
  if (i == 1){
    accuracy_df <- data.frame(marker_dense = marker_density[i], acc, G_corr = r)
  }else{
    accuracy_df <- rbind.data.frame(accuracy_df, data.frame(marker_dense = marker_density[i], acc, G_corr = r))
  }
}


# visualize accuracy
accuracy_df$text <- paste(accuracy_df$marker_dense, round(accuracy_df$G_corr,3), round(accuracy_df$PredCor,3), sep = '/')
ggplot(data=accuracy_df, aes(x=marker_dense, y=PredCor)) +
  geom_line() +
  geom_point() +
  xlab('Percentage of Markers Selected') +
  geom_text(aes(label = text),hjust=0, vjust=0)
  
