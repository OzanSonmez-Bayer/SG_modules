# Kinship Matrix

library(snpReady) # https://cran.r-project.org/web/packages/snpReady/vignettes/snpReady-vignette.html
library(tidyverse)

source('Code/Shared/genoFunctions.R')

# use pepper as an example for the calculation of Kinship Matrix

geno_dat <- read.delim(file = 'Pepper/Pepper_geno.txt')
map_dat <- read.delim(file = 'Pepper/Pepper_map.txt')

################################### Clean geno and map data ####################################

pos_chr <- map_dat %>% 
  select(X, CHROMOSOME, POSITION) %>% 
  mutate(X = as.character(X),
         CHROMOSOME = as.character(CHROMOSOME),
         POSITION = as.character(POSITION))


geno <- t(geno_dat)
colnames(geno) <- geno[1,]
X <- rownames(geno)
geno <- geno[-1,]

geno_filter <- filterfun(geno)
misind <- geno_filter[[1]][geno_filter[[1]] < 0.4] # remove individuals with more than 40% missing markers
mismrk <- geno_filter[[2]][geno_filter[[2]] < 0.25] # remove markers with more than 25% missing individuals. 

geno_c <- geno[, colnames(geno) %in% names(mismrk)]
geno_c[which(geno_c == "??", arr.ind = TRUE)] <- NA


# convert AA to 0, 1, 2
geno_convt <- gen2Additive(geno_c)
rownames(geno_convt) <- rownames(geno_c)
geno_filter <- filterfun(geno_convt)
misind <- geno_filter[[1]][geno_filter[[1]] < 0.6 ]
geno_cln <- geno_convt[rownames(geno_convt) %in% names(misind),] # 4589 individuals and 3715 markers


# Use snpReady to QC data
# start with geno raw data

geno[which(geno == "??", arr.ind = TRUE)] <- NA
geno_convt <- gen2Additive(geno) # convert nitrogenous to numeric
rownames(geno_convt) <- rownames(geno)

misind <- apply(geno_convt, 1, FUN = function(x) {sum(is.na(x))/length(x)})
mismark <- apply(geno_convt, 2, FUN = function(x) {sum(is.na(x))/length(x)})
hist(misind)
hist(mismark)
# Remove individuals with more than 60% missing values
# Remove markers with more than  75% missing obs
geno_ready <- raw.data(data = as.matrix(geno_convt), 
                       frame = 'wide', base = FALSE, 
                       sweep.sample = 0.7, 
                       call.rate = 0.25, 
                       maf = 0.1, 
                       imput = TRUE, 
                       imput.type = 'wright', 
                       outfile = '012')


# Wrigth method has 62% percent accuracy. Imputation basedon mean has intermediate accuracy as well. 

Mrk_cln <- geno_ready$M.clean


######################################## Compute G Matrix with snpReady ##################################

# Read in phenotypic data
#Use de-regressed GCA as input

# pheno <- read.csv('Pepper/BLUP.csv')
# gca_dat <- pheno %>% 
#   filter(mdl == 'gca' & trait == 'Fruit Shape') %>% 
#   mutate(pedigree = gsub('-', '.', pedigree))

gca_dat <- pred_df

# geno_sub <- geno_convet(which(rownames(geno_convet) %in% as.character(gca_dat$pedigree)))
rownames(Mrk_cln) <- gsub('\\.', '-', rownames(Mrk_cln))
Mrk_sub <- Mrk_cln[which(rownames(Mrk_cln) %in% as.character(gca_dat$Pedigree)),]
# remove monomorphic markers. These markers have the same genotype/form, which is not informative. 

mono_mrk <- apply(Mrk_sub, 2,  FUN = function(x){length(unique(x))})
mono_mrk <- names(mono_mrk)[which(mono_mrk == 1)]

Mrk_sub_1 <- Mrk_sub[,which(!colnames(Mrk_sub) %in% mono_mrk)]


# Compute G matrix with the select pedigrees
# Calculate G matrix
G <- G.matrix(M = Mrk_sub_1, method = 'VanRaden', format = 'long', plot = TRUE) # long format output is ready for Asreml
# Get additive matrix
G_mat <- G$Ga
# G matrix is usually not positive-definite
# Bending the matrix by diag(G) + 0.00001  or G_A = 0.99*G + 0.01*A

###################################### G-BLUP with Asreml ################################

library(asreml)
library(asremlPlus)

glup_afw <- asreml(fixed = dEBV ~ 1, 
                   random = ~Pedigree, 
                   ginverse = list(Pedigree = G_mat), 
                   data = gca_dat)

gpev <- predict(glup_afw, classify = 'pedigree')$predictions$pval


# Read in raw phenotypic data for parents
parents_pheno <- read.csv('Pepper/parental.csv')
parents_pheno$Pedigree <- gsub('-', '.', parents_pheno$Pedigree)

# generate a new kinship matrix

Mrk_sub_2 <- Mrk_cln[which(rownames(Mrk_cln) %in% as.character(parents_pheno$Pedigree)),]

mono_mrk_1 <- apply(Mrk_sub_2, 2,  FUN = function(x){length(unique(x))})
mono_mrk_1 <- names(mono_mrk_1)[which(mono_mrk_1 == 1)]

Mrk_sub_2 <- Mrk_sub_2[,which(!colnames(Mrk_sub_2) %in% mono_mrk_1)]

G_1 <- G.matrix(M = Mrk_sub_2, method = 'VanRaden', format = 'long', plot = FALSE)
G_mat_1 <- G_1$Ga

gblup_fruitshape <- asreml(fixed = FRSHP ~ 1, random = ~ giv(Pedigree, var = T), 
                           ginverse = list(Pedigree = G_mat_1), 
                           data = parents_pheno, 
                           na.method.Y = 'omit')

# GBLUP
pebv_fruitshape <- as.data.frame(predict(gblup_fruitshape, classify = 'Pedigree')$predictions$pval)

# Compare Raw with GBLUP

cor(pebv_fruitshape$predicted.value, parents_pheno$FRSHP[which(!is.na(parents_pheno$FRSHP))]) # 0.096

# May explore other package for GBLUP

######################################### rrBLUP ##############################################

library(rrBLUP)

load('Pepper/dEBV.Rdata')
# For rrBLUP, markers must be in {-1, 0, 1}
geno_rrblup <- geno_cln - 1

# Additive relationship matrix

geno_impute <- A.mat(geno_cln, impute.method = 'EM', return.imputed = TRUE)

Amatrix <- geno_impute$A # addtive relationship matrix

# Subset A matrix
Amatrix_sub <- Amatrix[which(rownames(Amatrix) %in% gca_dat$pedigree), which(colnames(Amatrix) %in% gca_dat$pedigree)]
dim(Amatrix_sub) # only 36 inbred lines

pred_df$Pedigree <- gsub('-', '.', pred_df$Pedigree)
Amatrix_sub <- Amatrix[which(rownames(Amatrix) %in% pred_df$Pedigree), which(colnames(Amatrix) %in% pred_df$Pedigree)]
dim(Amatrix_sub) # only 36 inbred lines

mdl_rrblup <- kin.blup(data = pred_df[pred_df$Pedigree %in% rownames(Amatrix_sub),], 
                       geno = 'Pedigree',
                       pheno = 'dEBV', 
                       K = Amatrix,
                       PEV = TRUE)
pred_ebv <- mdl_rrblup$pred #

cor(pred_df[pred_df$Pedigree %in% rownames(Amatrix_sub), "dEBV"], pred_ebv[names(pred_ebv) %in% rownames(Amatrix_sub)])  # correlation is -0.03

# This is an example: with genotype, we can predict phenotypic value. Amatrix consists non-phenotyped and phenotyped individuals. 

# Question: is it necessary to include all the SNPs in the matrix? 

######################################### Synbreed #################################

# install.packages('synbreed')
library(synbreed)

# Use synbreed package to generate synbreed object with genotypic and phenotyic data

# Match marker with Pedigree

# Step 1: Clean genotypic

geno_dat <- read.delim(file = 'Pepper/Pepper_geno.txt')
map_dat <- read.delim(file = 'Pepper/Pepper_map.txt')

pos_chr <- map_dat %>% 
  dplyr::select(X, CHROMOSOME, POSITION) %>% 
  dplyr::mutate(X = as.character(X),
         CHROMOSOME = as.character(CHROMOSOME),
         POSITION = as.character(POSITION))


geno <- t(geno_dat)
colnames(geno) <- geno[1,]
X <- rownames(geno)
geno <- geno[-1,]

geno[which(geno == "??", arr.ind = TRUE)] <- NA
geno_cln <- genClean(geno)
rownames(geno_cln) <- rownames(geno)

# Step 2: Phenotypic Data
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
  dplyr::select(X, CHROMOSOME, POSITION) %>% 
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

gblup_reml <- asreml(fixed = thk_dEBV ~ 1, random = ~giv(pedigree),
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


############################## Test GBLUP wtih New Lines #############################

# Mask phenotypic data, only use genotypic data

# Load crop ID
BaseInputBucketAddress <- "/shared/Veg_Phenotypes/Ref_Data/"
IDTableName <- paste0(BaseInputBucketAddress,'Pepper', "/",'Pepper',"_IDs.RData")
s3load(object = IDTableName, bucket = "genome-analytics-perm-space")

library(readxl)
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
  filter(as.character(pedigree) %in% rownames(geno_sub))

# Step 4: Create a Pedigree 
ped_dat <- pheno_sub %>% 
  select(pedigree) %>% 
  filter(!is.na(pedigree)) %>% 
  unique() %>% 
  mutate(ID = pedigree, 
         Par1 = 0, 
         Par2 = 0) %>%  
  select(ID, Par1, Par2)

# ped_dat[1,2:3] <- strsplit(as.character(ped_dat[1,1]), '\\+')[[1]] # hard coded line
ped_dat$gener <- rep(0, nrow(ped_dat))

ped_file <- create.pedigree(ped_dat$ID, ped_dat$Par1, ped_dat$Par2, ped_dat$gener)

plot(ped_file)


# Step 5: Create gpData object
gp <- create.gpData(geno = geno_sub, pheno = pheno_sub, pedigree = ped_file, map = map_sub)
summary(gp)

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
summary(Realized)
# 
# dimension                     82 x 82 
# rank                          81 
# range of off-diagonal values  -0.601138 -- 2.174691 
# mean off-diagonal values      -0.02352217 
# range of diagonal values      0.9531 -- 2.399 
# mean diagonal values          1.905295 

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
# range of off-diagonal values  -0.6081545 -- 2.197635 
# mean off-diagonal values      -0.02352759 
# range of diagonal values      0.9626 -- 2.382 
# mean diagonal values          1.905735 

plot(G, main = 'G Matrix')

par(mfrow = c(1,1))
boxplot(as.vector(G)~as.vector(Additive), xlab = 'A', ylab = 'G')

# Save G matrix in the format that is comptable with Asreml
Gtab <- write.relationshipMatrix(G, sorting = 'ASReml', type = 'ginv')
Atab <- write.relationshipMatrix(Additive, sorting = 'ASReml', type = 'ginv')


# Step 8: A-BLUP and G-BLUP
Ablup <- gpMod(gp_num, model = 'BLUP', kin = Additive, trait = 'LSM_thk')
summary(Ablup)

Gblup <- gpMod(gp_num, model = 'BLUP', kin = G, trait = 'LSM_thk')
summary(Gblup)
