# Genomic Selection: Elastic Net

library(aws.s3)
library(tidyverse)
library(caret)
library(glmnet)
library(jsonlite)
library(data.table)
library(Rcpp)

source('Code/Shared/misFunctions.R')
source('Code/Shared/ElasticNet.R')

AWS.Authorization('ycao1')

path <- 'ycao1/ITG_data/ProcTomato/'
bucketname <- 'genome-analytics-perm-space'

################ Download Genoptypic and Phenotypic Data ###############

# PGCA_PCM0_2014_2019_Deep.csv


pheno <- s3read_using(fread, 
                      object = 's3://genome-analytics-perm-space/ycao1/ITG_data/ProcTomato/BLUP/PGCA_PCM0_2014_2019_deep.csv')
# 4100 individuals
geno <- s3read_using(fread, 
                     object = 's3://genome-analytics-perm-space/ycao1/ITG_data/ProcTomato/Geno/imputedGeno.csv')
# 3302 individuals

############# Genotyped lines per Trait #############
trait_num <- c()
for (trt in unique(pheno$trait)){
  trait_num <- c(trait_num, length(intersect(geno$Pedigree, unique(pheno$PEDIGREE_NAME[pheno$trait == trt]))))
}

trait_df <- data.frame('Trait' = unique(pheno$trait), 
                        'Num of Genotyped Lines' = trait_num)
############# Calcuate DEBV ##############

pheno_debv <- pheno %>% 
  select(PEDIGREE_NAME, predicted.value, reliability, trait) %>% 
  mutate(dEBV = predicted.value/reliability) %>% 
  select(PEDIGREE_NAME, dEBV, trait)


############ Elastic Net with Random Search ###############

trtlist <- unique(pheno$trait)

accuracy <- data.frame()

outdir <- 'ProcTomato//GWS_Model/Imputed/ENET_random/'

for (trtname in trtlist){
  rsquare_df <- TrainModels_randomsearch(pheno = pheno_debv, 
                                         geno = geno,
                                         trt = trtname,
                                         outdir = outdir)
  accuracy <- rbind(accuracy, rsquare_df)
}

####################### Optimze lambda Give Alpha #############
set.seed(20200416)

outdir <- 'ProcTomato/GWS_Model/Imputed/ENET/'
alpha_list <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0)
for (trtname in trtlist){
  print(trtname)
  TrainModels_predefined(pheno = pheno_debv, geno = geno, trt = trtname, outdir = outdir, alpha_list = alpha_list)
}


############################################### QA/QC Marker Data ##################################

library(snpReady)

geno_matrix <- as.matrix(geno[,-c('ProgenyGermID', 'Pedigree')])
rownames(geno_matrix) <- geno$Pedigree

# Remove markers with MAF less than 0.05

geno_cln <- raw.data(geno_matrix, frame = 'wide', hapmap = NULL, base = FALSE, outfile = '012', imput = FALSE)
geno_ready <- geno_cln$M.clean
geno_report <- geno_cln$report


############## Elastic Net with Random Search and QCed data #######################

outdir <- 'ProcTomato//GWS_Model/Imputed/QA/ENET_random/'

geno_ready <- data.frame(geno[, c('ProgenyGermID', 'Pedigree')], geno_ready)

accuracy <- data.frame()

for (trtname in trtlist){
  print(trtname)
  rsquare_df <- TrainModels_randomsearch(pheno = pheno_debv, 
                                         geno = geno_ready,
                                         trt = trtname,
                                         outdir = outdir)
  accuracy <- rbind(accuracy, rsquare_df)
}


outdir <- 'ProcTomato//GWS_Model/Imputed/QA/ENET_predefined/'
alpha_list <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001, 0)
for (trtname in trtlist){
  print(trtname)
  TrainModels_predefined(pheno = pheno_debv, geno = geno_ready, trt = trtname, outdir = outdir, alpha_list = alpha_list)
}

