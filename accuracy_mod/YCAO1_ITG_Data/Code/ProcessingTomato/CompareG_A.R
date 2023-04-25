# Compare A matrix with G matrix that is calculated wtih imputed genotype

library(tidyverse)
library(aws.s3)
library(ggplot2)
library(jsonlite)
library(data.table)
library(bit64)
library(AGHmatrix)
library(dplyr)

source('Code/Shared/misFunctions.R')

AWS.Authorization('ycao1')

path <- 'ycao1/ITG_data/ProcTomato/'
bucketname <- 'genome-analytics-perm-space'

A_2gen <- s3read_using(fread, object = paste0(path, 'Amat/Additive_2gen.csv'), bucket = bucketname)
A_5gen <- s3read_using(fread, object = paste0(path, 'Amat/Additive_5gen.csv'), bucket = bucketname)
  
geno <- s3read_using(fread, object = paste0(path, 'Geno/imputedGeno.csv'), bucket = bucketname)

rownames(A_2gen) <- colnames(A_2gen)
A_2gen <- as.matrix(A_2gen)
# rownames(A_5gen) <- colnames(A_5gen) # Some issues with A 5gen pedigrees

# G matrix calculation

geno_ready <- geno %>% 
  dplyr::select(-c('ProgenyGermID', 'Pedigree'))
rownames(geno_ready) <- geno$Pedigree
geno_ready <- as.matrix(geno_ready)

G_mat <- Gmatrix(geno_ready, method = 'VanRaden')


# Match individuals in pheno and geno

ped_list <- intersect(rownames(A_2gen), rownames(G_mat))
length(ped_list)

G_sub <- G_mat[which(rownames(G_mat) %in% ped_list), which(colnames(G_mat) %in% ped_list)]
A_sub <- A_2gen[which(rownames(A_2gen) %in% ped_list), which(colnames(A_2gen) %in% ped_list)]

G_sub <- G_sub[rownames(A_sub), colnames(A_sub)]


hist(G_sub[diag(G_sub)])
range(G_sub[diag(G_sub)])

hist(A_sub[diag(A_sub)])
range(A_sub[diag(A_sub)])
