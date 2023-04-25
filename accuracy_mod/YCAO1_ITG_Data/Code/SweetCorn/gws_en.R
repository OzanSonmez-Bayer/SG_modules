# Sweet Corn Elastic Net
library(aws.s3)
library(tidyverse)
library(caret)
library(glmnet)
library(jsonlite)
library(data.table)
library(Rcpp)

if (length(grep('DOMINO', names(Sys.getenv()))) > 0){
  setwd('/repos/ITG_Data')
}


source('Code/Shared/ElasticNet.R')
source('Code/Shared/misFunctions.R')

AWS.Authorization('ycao1')

path <- 'ycao1/ITG_data/SweetCorn/'
bucketname <- 'genome-analytics-perm-space'


pheno <- s3read_using(fread, 
                      object = paste(path, 'BLUP/PGCA_Proc_2014_2018_deep.csv', sep = ''), 
                      bucket = bucketname)
# 8396 individuals
geno <- s3read_using(fread, 
                     object = paste(path, 'Geno/imputedGeno.csv', sep = ''),
                     bucket = bucketname)
# 8962 individuals


############# Genotyped lines per Trait #############
trait_num <- c()
for (trt in unique(pheno$trait)){
  trait_num <- c(trait_num, length(intersect(geno$Pedigree, unique(pheno$PEDIGREE_NAME[pheno$trait == trt]))))
}

trait_df <- data.frame('Trait' = unique(pheno$trait), 
                       'Num of Genotyped Lines' = trait_num)
trait_df

############# Calcuate DEBV ##############

pheno_debv <- pheno %>% 
  select(PEDIGREE_NAME, predicted.value, reliability, trait) %>% 
  mutate(dEBV = predicted.value/reliability) %>% 
  select(PEDIGREE_NAME, dEBV, trait)

rm(pheno)

############ Elastic Net ################

# dir.create(path = 'Corn/GWS_Model/ENET_random')

if (length(grep('DOMINO', names(Sys.getenv()))) > 0){
  outdir <- '/mnt/GWS/Corn/ENET_random/'
}else{
  outdir <- '/Corn/GWS_Model/ENET_random/'
}


trtlist <- unique(pheno_debv$trait)

accuracy <- data.frame()

for (trtname in trtlist){
  print(trtname)
  rsquare_df <- TrainModels_randomsearch(pheno = pheno_debv, 
                                         geno = geno,
                                         trt = trtname,
                                         outdir = outdir)
  accuracy <- rbind(accuracy, rsquare_df)
}


outdir <- '/mnt/GWS/Corn/Processing/ENET_predefined/'
alpha_list <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0)
for (trtname in trtlist){
  print(trtname)
  TrainModels_predefined(pheno = pheno_debv, geno = geno, trt = trtname, outdir = outdir, alpha_list = alpha_list)
}


