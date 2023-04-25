install.packages('/mnt/packages/glmnet_2.0-18.tar.gz', type = 'source', repos = NULL)

library(aws.s3)
library(tidyverse)
library(caret)
library(glmnet)
library(jsonlite)
library(data.table)
library(Rcpp)
library(httr)


source('vaultCredentials.R')
source('/repos/BLUPF90/BLUPF90/funcStore_deep_ITG.R')
source('/repos/ITG_Data/Code/Shared/ElasticNet.R')

loadcrential()

args = commandArgs(trailingOnly = TRUE)

market <- args[1]
outputfile <- args[2]

if (market == 'fresh'){
  alldata <- s3readRDS(object = 'ycao1/ITG/fresh_data.rds',
                       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
  trait <- c('HSC', 'SC_TF', 'SDV', 'ELEN', 'QUAL', 'SC_PR', 'EDIA')
  
}else if(market == 'process'){
  alldata <- s3readRDS(object = 'ycao1/ITG/process_data.rds',
                       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
  trait <- c('HSC', 'RTLR', 'SDV', 'SC_HE', 'S50D')
}else{
  alldata <- s3readRDS(object = 'ycao1/ITG/pheno.rds',
                       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
  trait <- c('SC_TF', 'EDIA', 'ELEN', 'HSC', 'PA', 'SC_HE', 'SC_SE', 'SC_PR', 'SDV', 'ESHK', 'TILLR', 'RTLR', 'QUAL',
             'SC_FG', 'SC_FR', 'S50', 'S50_BE', 'LDSR_SETOTU', 'LDSR_PUCCSO')
}



geno <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds',
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

# load('/domino/datasets/local/hybridGeno/g_hybrid.RData')



alpha_list <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0)

for (trt in trait){
  print(trt)
  trainingModel(pheno = alldata, 
                geno = geno, 
                trait = trt, 
                alpha_list = alpha_list, 
                outputfile = outputfile, 
                validation_only = F)
}




# ############ Check Results ###############
# # # 
# load('/mnt/Corn/ElasticNet/testing/SC_FR.rds')
# results_df
