install.packages('/mnt/packages/glmnet_2.0-18.tar.gz', type = 'source', repos = NULL)

library(aws.s3)
library(tidyverse)
library(caret)
library(glmnet)
library(jsonlite)
library(data.table)
library(Rcpp)
library(httr)


# source('vaultCredentials.R')
source('/repos/BLUPF90/BLUPF90/funcStore_deep_ITG.R')
source('/repos/ITG_Data/Code/Shared/ElasticNet.R')

loadcrential()

args = commandArgs(trailingOnly = TRUE)

# market <- args[1]
outputfile <- args[1]
validation <- args[2]
testing_season <- args[3]


# Get Training data

######################### Get Pheno Training #######################

pheno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/pheno_training.rds',
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


######################### Get Hybrid Geno Training #######################


# geno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/hybridGeno_training.rds',
                   # bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')



geno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/geno_all.rds',
                  bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

######################## Model Training #################################
trait <- c('FBRIX', 'AVGHBRIX', 'AVJB', 'FQUAL', 'FRFRMH', 'FZUNI', 'LBRIX', 'MAT', 'MATR', 'OST', 'PLCMP', 'PLTVG', 'PYLDPA', 'SAMPUW')

alpha_list <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0)

# outputfile <- 'ProcTomato/ElasticNet/'

for (trt in trait){
  print(trt)
  trainingModel(pheno = pheno, 
                geno = geno, 
                trait = trt, 
                alpha_list = alpha_list, 
                outputfile = outputfile, 
                validation_only = validation, 
                testing_season = testing_season)
}

# 
# tst %>% 
#   group_by(OBSRVTN_REF_CD) %>% 
#   summarise(n = length(unique(PEDIGREE_NAME)))
