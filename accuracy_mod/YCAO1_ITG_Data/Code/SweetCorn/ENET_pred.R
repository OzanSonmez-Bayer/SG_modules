install.packages('/mnt/packages/glmnet_2.0-18.tar.gz', type = 'source', repos = NULL)

library(aws.s3)
library(tidyverse)
library(caret)
library(glmnet)
library(jsonlite)
library(data.table)
library(Rcpp)
library(httr)
library(parallel)
library(doParallel)
registerDoParallel(16)

source('vaultCredentials.R')
source('/repos/BLUPF90/BLUPF90/funcStore_deep_ITG.R')
source('/repos/ITG_Data/Code/Shared/ElasticNet.R')

loadcrential()


geno <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds',
                  bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

alldata <- s3readRDS(object = 'ycao1/ITG/Corn/pheno.rds',
                     bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

s3load(object = 'Corn/Corn_IDs.RData', bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")
CropIDs$M.GERMPLASM.PEDIGREE <- as.character(CropIDs$M.GERMPLASM.PEDIGREE)

s3load(object = 'GPC_Data/Corn/Corn_GPC.RData', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data')
GCP_TableView_noProg$Pedigree <- as.character(GCP_TableView_noProg$Pedigree)


# trait <- 'EDIA'
testing_season <- '2020:04'
# alpha_val <- 0
trait_alpha <- data.frame(trait = c('SC_TF', 'EDIA', 'ELEN', 'HSC', 'PA', 'SC_HE', 'SC_SE', 'SC_PR', 'SDV', 'ESHK',
                                    'TILLR', 'RTLR', 'QUAL', 'SC_FG', 'SC_FR', 'S50', 'S50_BE', 'LDSR_SETOTU', 'LDSR_PUCCSO'), 
                          alpha = c(0.05, 0, 0.001, 0, 1, 0, 0.01, 0, 1, 0, 0.4, 0, 1, 0.6, 0, 0, 0.001, 0, 0))
hybrid_ls <- unique(as.character(alldata$PEDIGREE_NAME))

for (trait in trait_alpha$trait){
  print(trait)
  
  alpha_val <- trait_alpha$alpha[trait_alpha$trait == trait]
  
  pheno <- alldata %>% 
    filter(OBSRVTN_REF_CD == trait & GROWSEASON != testing_season) %>% 
    select(PEDIGREE_NAME, TRAIT_VALUE) %>% 
    mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
    group_by(PEDIGREE_NAME) %>% 
    summarize(lsmean = mean(TRAIT_VALUE, na.rm = TRUE))
  
  pheno_geno <- pheno %>% 
    inner_join(geno, by = c('PEDIGREE_NAME' = 'Pedigree'))
  
  pheno_w_testingseason <- alldata %>% 
    filter(OBSRVTN_REF_CD == trait) %>% 
    select(PEDIGREE_NAME, TRAIT_VALUE) %>% 
    mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
    group_by(PEDIGREE_NAME) %>% 
    summarize(lsmean = mean(TRAIT_VALUE, na.rm = TRUE))
  
  pheno_geno_w_ts <- pheno_w_testingseason %>% 
    inner_join(geno, by = c('PEDIGREE_NAME' = 'Pedigree'))
  
  set.seed(2020)
  
  training_indx <- createDataPartition(pheno_geno$lsmean, times = 1, p = 0.9, list = FALSE)
  
  training_set <- pheno_geno[training_indx, ]
  testing_set <- pheno_geno[-training_indx, ]
  
  
  y_train <- training_set %>% 
    select(lsmean) %>% 
    as.matrix()
  X_train <- training_set %>% 
    select(-c(lsmean, PEDIGREE_NAME, GermID)) %>% 
    as.matrix()
  
  y_test <- testing_set %>% 
    select(lsmean) %>% 
    as.matrix()
  X_test <- testing_set %>% 
    select(-c(lsmean, PEDIGREE_NAME, GermID)) %>% 
    as.matrix()
  
  mdl_tuning <- cv.glmnet(X_train, y_train, nfolds = 5, alpha = alpha_val,parallel = TRUE)
  mdl <- glmnet(X_train, y_train, lambda = mdl_tuning$lambda.min, alpha = alpha_val)
  
  
  X_hybrid <- pheno_geno_w_ts %>%
    select(-c(lsmean, PEDIGREE_NAME, GermID)) %>%
    as.matrix()

  y_hybrid <- pheno_geno_w_ts %>%
    select(lsmean) %>%
    as.matrix()

  y_hat_hybrid <- mdl %>% predict(X_hybrid) %>% as.vector()
  rsq <- cor(y_hybrid, y_hat_hybrid)
  rsq

  GEBV_hybrid <- data.frame(Pedigree = pheno_geno_w_ts$PEDIGREE_NAME, 
                            GEBV = y_hat_hybrid, 
                            Raw_LSMean = pheno_geno_w_ts$lsmean)

  # Add Origin and HBC

  GEBV_hybrid <- GEBV_hybrid %>%
    left_join(CropIDs[, c('M.GERMPLASM.PEDIGREE', 'M.GERMPLASM.ORIGIN')], by = c('Pedigree' = 'M.GERMPLASM.PEDIGREE')) %>%
    rename(origin = M.GERMPLASM.ORIGIN) %>%
    unique()  %>%
    left_join(GCP_TableView_noProg[, c('Pedigree', 'HBC')], by = 'Pedigree') %>%
    unique()

  s3write_using(GEBV_hybrid, FUN = write.csv,
                object = paste0('ycao1/ITG/Corn/ENET/', trait, '_GEBV_hybrid.csv'),
                bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

  rm(GEBV_hybrid)
  rm(X_hybrid)
  rm(y_hybrid)
  gc()

  # # inbred_ls <- unique(c(as.character(alldata[alldata$OBSRVTN_REF_CD == trait, 'P1']),
  # #                       as.character(alldata[alldata$OBSRVTN_REF_CD == trait, 'P2'])))
  # 
  # # X_inbred <- geno %>% 
  # #   filter(as.character(Pedigree) %in% inbred_ls) %>% 
  # #   select(-c(Pedigree, GermID)) %>% 
  # #   as.matrix()
  
  # X_inbred <- geno %>% 
  #   filter(!as.character(Pedigree) %in% hybrid_ls) %>% 
  #   select(-c(Pedigree, GermID)) %>% 
  #   as.matrix()
  # 
  # dim(X_inbred)
  # 
  # y_inbred <- mdl %>% 
  #   predict(X_inbred) %>% 
  #   as.vector()
  # # GEBV_inbred <- data.frame(Pedigree = geno$Pedigree[which(geno$Pedigree %in% inbred_ls)], GEBV = y_inbred)
  # 
  # GEBV_inbred <- data.frame(Pedigree = geno$Pedigree[which(!geno$Pedigree %in% hybrid_ls)], GEBV = y_inbred)
  # 
  # 
  # GEBV_inbred <- GEBV_inbred %>% 
  #   left_join(CropIDs[, c('M.GERMPLASM.PEDIGREE', 'M.GERMPLASM.ORIGIN')], by = c('Pedigree' = 'M.GERMPLASM.PEDIGREE')) %>% 
  #   rename(origin = M.GERMPLASM.ORIGIN) %>% 
  #   unique()  %>%
  #   left_join(GCP_TableView_noProg[, c('Pedigree', 'HBC')], by = 'Pedigree') %>%
  #   unique()
  # 
  # s3write_using(GEBV_inbred, FUN = write.csv, 
  #               object = paste0('ycao1/ITG/Corn/ENET/', trait, '_GEBV_inbred_all.csv'), 
  #               bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
  # 
  # rm(GEBV_inbred)
  # gc()
}


### Collect Results ####

GEBV_hybrid <- data.frame()

GEBV_inbred <- data.frame()

for (tt in trait_alpha$trait){
  
  # Hybrid
  temp_1 <- s3read_using(fread,
               object = paste0('ycao1/ITG/Corn/ENET/', tt, '_GEBV_hybrid.csv'),
               bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
  temp_1 <- temp_1[,-1]

  colnames(temp_1)[which(colnames(temp_1) %in% c('GEBV', 'Raw_LSMean'))] <- paste0(tt, '_', colnames(temp_1)[which(colnames(temp_1) %in% c('GEBV', 'Raw_LSMean'))])

  if(nrow(GEBV_hybrid) == 0){
    GEBV_hybrid <- temp_1
  }else{
    GEBV_hybrid <- GEBV_hybrid %>%
      full_join(temp_1, by = c('Pedigree', 'origin', 'HBC'))
  }
  
  # Inbred
  
  # temp_2 <- s3read_using(fread, 
  #                        object = paste0('ycao1/ITG/Corn/ENET/', tt, '_GEBV_inbred_all.csv'), 
  #                        bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
  # temp_2 <- temp_2[,-1]
  # 
  # colnames(temp_2)[which(colnames(temp_2) %in% c('GEBV'))] <- paste0(tt, '_', colnames(temp_2)[which(colnames(temp_2) %in% c('GEBV'))])
  # 
  # if(nrow(GEBV_inbred) == 0){
  #   GEBV_inbred <- temp_2
  # }else{
  #   GEBV_inbred <- GEBV_inbred %>% 
  #     full_join(temp_2, by = c('Pedigree', 'origin', 'HBC'))
  # }
  # 
  rm(temp_1)
  # rm(temp_2)
  gc()
  
}


s3write_using(GEBV_hybrid, FUN = write.csv, 
              object = 'ycao1/ITG/Corn/ENET/GEBV_hybrid_alltraits.csv',
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
s3write_using(GEBV_inbred, FUN = write.csv, 
              object = 'ycao1/ITG/Corn/ENET/GEBV_inbred_alltraits.csv',
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
s3write_using(GEBV_inbred, FUN = write.csv, 
              object = 'ycao1/ITG/Corn/ENET/GEBV_inbred_alltraits_allinbreds.csv',
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')