# Compare A matrix with G matrix that is calculated wtih imputed genotype

library(tidyverse)
library(aws.s3)
library(ggplot2)
library(jsonlite)
library(data.table)
library(bit64)
install.packages('/mnt/packages/AGHmatrix_0.0.4.tar.gz', type = 'source', repos = NULL)

library(AGHmatrix)
library(dplyr)
library(httr)
library(asreml)
library(asremlPlus)

# source('Code/Shared/misFunctions.R')

source('Credential/vaultCredentials.R')
source('/repos/BLUPF90/BLUPF90/funcStore_deep_ITG.R')
source('/repos/ITG_Data/Code/Shared/ElasticNet.R')

loadcrential()

Crop <- 'Corn'
source('/repos/ITG_Data/Code/Shared/AmatCalc.R')




geno <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds',
                  bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
alldata <- s3readRDS(object = 'ycao1/ITG/Corn/pheno.rds',
                     bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

# rownames(A_2gen) <- colnames(A_2gen)
# A_2gen <- as.matrix(A_2gen)
# rownames(A_5gen) <- colnames(A_5gen) # Some issues with A 5gen pedigrees

# G matrix calculation
geno <- geno[!duplicated(geno$Pedigree), ]
geno_ready <- geno %>% 
  dplyr::select(-c('GermID', 'Pedigree'))
rownames(geno_ready) <- geno$Pedigree
geno_ready <- as.matrix(geno_ready)

G_mat <- Gmatrix(geno_ready, method = 'VanRaden')


# Calculate A matrix 

data_sub <- alldata %>% 
  filter(OBSRVTN_REF_CD == 'ELEN')

A_5gen <- Amat_5gen(data_sub, 5)

############# Compare G and A ################

# Match individuals in pheno and geno

A_5gen <- s3readRDS(object = 'ycao1/ITG/Corn/A_5gen.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
G_mat <- s3readRDS(object = 'ycao1/ITG/Corn/GRM.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

ped_list <- intersect(rownames(A_5gen), rownames(G_mat))
length(ped_list)

G_sub <- G_mat[which(rownames(G_mat) %in% ped_list), which(colnames(G_mat) %in% ped_list)]
A_sub <- A_5gen[which(rownames(A_5gen) %in% ped_list), which(colnames(A_5gen) %in% ped_list)]

G_sub <- G_sub[rownames(A_sub), colnames(A_sub)]


hist(G_sub[diag(G_sub)])
range(G_sub[diag(G_sub)])

hist(A_sub[diag(A_sub)])
range(A_sub[diag(A_sub)])
