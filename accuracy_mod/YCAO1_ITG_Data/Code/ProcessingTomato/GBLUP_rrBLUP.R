# GBLUP for Processing Tomato

library(readxl)
library(tidyverse)
library(snpReady)
library(synbreed)
library(ggplot2)
library(asreml)
library(asremlPlus)
library(aws.s3)
library(jsonlite)
library(Matrix)

source('Code/Shared/BLUPFunctions.R')
source('Code/Shared/genoFunctions.R')
source('Code/Shared/misFunctions.R')
source('Code/Shared/CGS.R')

################# Import Raw Taqman Data #################

procTom_taqman <- read.table('ProcTomato/9Z_TaqManFromCGS_FilteredMarkers.tab', 
                             sep = '\t',
                             fill = TRUE,
                             header = TRUE)
procTom_taqman <- procTom_taqman %>% 
  filter(!is.na(Pedigree))

dim(procTom_taqman)
miss_mrk <- apply(procTom_taqman[,-c(1:2)], MARGIN = 1, FUN = function(x){sum(is.na(x))})
miss_ind <- apply(procTom_taqman[,-c(1:2)], MARGIN = 2, FUN = function(x){sum(is.na(x))})

# remove individuals with 285 missing markers
procTom_taqman_1 <- procTom_taqman[-which(miss_mrk == 287),] # 2607 individuals left
miss_ind <- apply(procTom_taqman_1[,-c(1:2)], MARGIN = 2, FUN = function(x){sum(is.na(x))})
procTom_taqman_2 <- procTom_taqman_1[, -which(miss_ind == 2599)]
geno <- procTom_taqman_2[,-1]


procTom_taqman_2[1:5,1:5]

taqman_full <- procTom_taqman_2 %>% 
  dplyr::select(-Inventory_SampleID) %>% 
  unique()
rownames(taqman_full) <- taqman_full$Pedigree
taqman_full <- taqman_full %>% 
  dplyr::select(-Pedigree)

taqman_full_1 <- gen2Additive(taqman_full)
rownames(taqman_full_1) <- rownames(taqman_full)

G_rrBLUP <- rrBLUP::A.mat((taqman_full_1 - 1), 
                          min.MAF = 0.01,  
                          impute.method = 'mean', 
                          tol = 0.02, 
                          return.imputed = TRUE)

G_imputed <- G_rrBLUP$imputed
hist(diag(G_rrBLUP$A))
range(diag(G_rrBLUP$A)) 



############### Import GCA from Deep Pedigree #############

AWS.Authorization('ycao1')
path <- 'ycao1/ITG_data/ProcTomato/'
bucketname <- 'genome-analytics-perm-space'
infile <- paste0('s3://', bucketname, '/', path, 'BLUP/', 'PGCA_PCM0_2014_2019_deep.csv')
csvfile <- get_object(infile)
csvfile_1 <- rawToChar(csvfile)
con <- textConnection(csvfile_1)
pheno <- read.csv(con)
close(con)


# taqman_full_sub <- taqman_full[which(rownames(taqman_full) %in% as.character(unique(pheno$PEDIGREE_NAME))),]
# taqman_full_sub_012 <- gen2Additive(taqman_full_sub)
# rownames(taqman_full_sub_012) <- rownames(taqman_full_sub)
# 
# 
# G_rrBLUP <- rrBLUP::A.mat((taqman_full_sub_012 - 1), 
#                           min.MAF = 0.01, 
#                           impute.method = 'mean', 
#                           tol = 0.02)
# hist(diag(G_rrBLUP))
# range(diag(G_rrBLUP)) 



############## Match Geno with Pheno #########################

# FBRIX
# MAT, FRFRMH, PLTVG, PLCMP, FQUAL, PYLDPA, SAMPUW, OST, AVJB, LBRIX, FBRIX, AVGHBRIX, FZUNI, MATR, EARLY
trt <- 'SAMPUW'
pheno_df <- pheno %>% 
  filter(trait == trt) %>% 
  mutate(dEBV = predicted.value/reliability)

length(intersect(as.character(pheno_df$PEDIGREE_NAME), rownames(G_rrBLUP))) # 496 inbreds have been genotyped

pheno_sub <- pheno_df %>% 
  filter(as.character(PEDIGREE_NAME) %in% rownames(G_rrBLUP)) %>%
  dplyr::select(PEDIGREE_NAME, dEBV)
pheno_sub <- pheno_sub[!duplicated(pheno_sub$PEDIGREE_NAME),]
rownames(pheno_sub) <- pheno_sub$PEDIGREE_NAME
pheno_sub$PEDIGREE_NAME <- as.character(pheno_sub$PEDIGREE_NAME)
# pheno_sub <- pheno_sub %>% 
#   dplyr::select(-PEDIGREE_NAME)


G_sub <- G_rrBLUP[rownames(G_rrBLUP) %in% pheno_sub$PEDIGREE_NAME, colnames(G_rrBLUP) %in% pheno_sub$PEDIGREE_NAME]

realizedPD <- nearPD(G_sub, keepDiag = T)
G_base <- matrix(realizedPD[[1]]@x, nrow = realizedPD[[1]]@Dim[1])
G_base <- G_base + diag(0.01, nrow = nrow(G_base))
attr(G_base, 'dimnames') <- realizedPD[[1]]@Dimnames
class(G_base) <- 'relationshipMatrix'
str(G_base)
summary(G_base)

G_base <- G_base[pheno_sub$PEDIGREE_NAME, pheno_sub$PEDIGREE_NAME]


plot(G_base)



G_inv <- write.relationshipMatrix(G_base, sorting = 'ASReml', type = 'ginv')
head(attr(G_inv, 'rowNames'))
names(G_inv) <- c('row', 'column', 'coefficient')
head(G_inv)


############### GBLUP ###############

GBLUP_asreml <- asreml(fixed = dEBV ~ 1,
                       random = ~giv(PEDIGREE_NAME, var = T), 
                       ginverse = list(PEDIGREE_NAME = G_inv), 
                       rcov = ~units, 
                       data = pheno_sub, 
                       na.method.X = 'include',
                       verbose = TRUE,
                       workspace = 320e+6, pworkspace = 320e+6)
pred_asreml <- as.data.frame(predict(GBLUP_asreml, classify = 'PEDIGREE_NAME')$predictions$pval)
cor(pheno_sub$dEBV, pred_asreml$predicted.value)

varcomp_tab <- summary(GBLUP_asreml)$varcomp
varcomp_tab$component[1]/sum(varcomp_tab$component)



# install.packages('sommer')
library(sommer)
ans <- mmer(dEBV~1,
            random=~vs(PEDIGREE_NAME,Gu=G_rrBLUP),
            rcov=~units,
            data=pheno_sub)

pin(ans, h2 ~ (V1) / ( V1+V2) )


################### ElasticNet for GBLUP #######################

library(glmnet)
library(caret)

trtlist <- unique(pheno$trait)
accuracy_list <- list()
coef_list <- list()
for (trt in trtlist){
  # trt <- 'EARLY'
  print(trt)
  pheno_df <- pheno %>% 
    filter(trait == trt) %>% 
    mutate(dEBV = predicted.value/reliability) %>% 
    dplyr::select(PEDIGREE_NAME, dEBV)
  
  # Match Pheno with Geno
  
  pheno_geno <- as.data.frame(G_imputed) %>% 
    mutate(PEDIGREE_NAME = rownames(.)) %>% 
    inner_join(pheno_df,., by = 'PEDIGREE_NAME')
  
  pheno_geno <- pheno_geno[!duplicated(pheno_geno$PEDIGREE_NAME), ]
  
  
  rownames(pheno_geno) <- pheno_geno$PEDIGREE_NAME
  
  pheno_geno_1 <- pheno_geno %>% 
    dplyr::select(-PEDIGREE_NAME)
  
  
  # Train model 
  
  set.seed(20200416)
  
  # inTraining <- createDataPartition(pheno_geno_1$dEBV, p = .8, list = FALSE)
  # training <- pheno_geno_1[ inTraining,]
  # testing  <- pheno_geno_1[-inTraining,]
  # 
  # x_train <- as.matrix(training[,-1])
  # y_train <- training$dEBV
  # 
  # x_test <- as.matrix(testing[,-1])
  # y_test <- testing$dEBV
  
  
  x_train <- as.matrix(pheno_geno_1[,-1])
  y_train <- as.matrix(pheno_geno_1$dEBV)
  fit.lasso <- glmnet(x_train, y_train, alpha=1)
  fit.ridge <- glmnet(x_train, y_train, alpha=0)
  fit.elnet <- glmnet(x_train, y_train, alpha=.5)
  
  for (i in 0:10) {
    assign(paste("fit", i, sep=""), cv.glmnet(x_train, y_train, type.measure="mse", 
                                              alpha=i/10))
  }
  
  plot(fit.lasso, xvar="lambda")
  plot(fit10, main="LASSO")

  plot(fit.ridge, xvar = 'lambda')
  plot(fit0, main = 'Ridge')

  plot(fit.elnet, xvar = 'lambda')
  plot(fit5, main = 'Elastic Net')
  
  
  predicted_gt <- data.frame(dEBV = y_train)
  
  for (i in 0:10){
    temp <- get(paste('fit', i, sep = ''))
    # assign(paste('pred', i, sep = '_'), predict(temp, s = temp$lambda.1se, newx = x_test))
    temp_mdl <- glmnet(x_train, y_train, lambda = temp$lambda.min)
    # temp_pred <- temp_mdl %>% predict(x_train) %>% as.vector()
    assign(paste(trt, 'coeff', i, sep = '_'), coef.glmnet(temp_mdl)@x)
    coef_list[[paste(trt, i, sep = '_')]] <- get(paste(trt, 'coeff', i, sep = '_'))
    temp_pred <- predict(temp, s = temp$lambda.min, newx = x_train)
    predicted_gt <- data.frame(predicted_gt, temp_pred)
    
  }
  
  colnames(predicted_gt) <- c('dEBV', paste('alpha', c(0:10)/10, sep = '_'))
  
  # Calculate MSE
  mse_df <- colMeans((predicted_gt[, -1] - predicted_gt[,1]) * (predicted_gt[,-1] - predicted_gt[, 1]))
  
  RMSE_val <- apply(predicted_gt[, -1], MARGIN = 2, FUN = function(x){RMSE(x, predicted_gt[, 1])})
  R2_val <- apply(predicted_gt[, -1], MARGIN = 2, FUN = function(x){cor(x, predicted_gt[, 1])})
  
  
  accuracy_df <- data.frame(RMSE_val = RMSE_val, correlation = R2_val)
  accuracy_df
  
  accuracy_list[[trt]] <- accuracy_df
  
}


str(coef.glmnet(glmnet(x_train, y_train, lambda = fit10$lambda.min)))


# Can we use raw data as y and add grow seasons, trep, and field as covariate? 
