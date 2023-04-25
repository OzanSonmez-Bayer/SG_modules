library(readxl)
library(tidyverse)
library(aws.s3)
library(plyr)
library(dplyr)

library(asreml)
library(asremlPlus)
library(doBy)
library(jsonlite)

source('Code/Shared/misFunctions.R')
source('Code/Shared/BLUPFunctions.R')
AWS.Authorization('ycao1')

Crop <- 'Corn'

fresh_data <- data.frame()
proc_data <- data.frame()

path <- 'ycao1/ITG_data/SweetCorn/Data/'
bucket <- 'genome-analytics-perm-space'

############## Pull Clean Raw Data from S3 #################

file_list <- get_bucket(bucket = bucket, prefix = paste(path, 'Raw/', sep = ""))
file_list <- file_list[grepl('_Raw', file_list)]


for (i in 1:length(file_list)){
  print(file_list[[i]]$Key)
  infile <- paste0('s3://','genome-analytics-perm-space/', file_list[[i]]$Key)
  csvfile <- get_object(infile)
  csvfile_1 <- rawToChar(csvfile)
  con <- textConnection(csvfile_1)
  temp_dat <- read.csv(con)
  close(con)
  if (grepl('_Fresh_', file_list[[i]]$Key)){
    fresh_data <- plyr::rbind.fill(fresh_data, temp_dat)
  }else{
    proc_data <- plyr::rbind.fill(proc_data, temp_dat)
  }
}

###################### Generate 2-gen A matrix #####################

Amat_2gen <- function(dat) {
  
  if (!'P1' %in% colnames(dat)){
    print("parentage pedigree doesn't exit")
    dat <- dat %>% 
      separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
    p1 <- unlist(lapply(dat$P1, FUN = pedSingle))
    p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    dat$P1 <- p1
    dat$P2[dat$P2=='']=NA
    p2 <- unlist(lapply(dat$P2, FUN = pedSingle))
    p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    dat$P2 <- p2
  }
  raw_dat_sub <- dat %>% 
    mutate(P1 = as.character(P1), 
           P2 = as.character(P2),
           trait = OBSRVTN_REF_CD, 
           TRAIT_VALUE = as.numeric(TRAIT_VALUE), 
           GROWSEASON = as.character(GROWSEASON),
           pedigree = as.character(PEDIGREE_NAME))
  raw_dat_sub[is.na(raw_dat_sub$P1), 'P2'] <- NA
  raw_dat_sub[is.na(raw_dat_sub$P2), 'P1'] <- NA
  
  
  ped_file <- raw_dat_sub %>% 
    dplyr::select(PEDIGREE_NAME, P1, P2) %>% 
    distinct_all()
  
  # Two generation A matrix
  
  ped_file$Selfing <- rep(0, rep = nrow(ped_file))
  A_2gen <- asreml.Ainverse(ped_file, fgen = c('Selfing',5))
  A_inv_2gen <- A_2gen$ginv
  A_inv_sparse_2gen <- asreml.sparse2mat(A_inv_2gen)
  Additive_2gen <- round(solve(A_inv_sparse_2gen),2) 
  
  
  colnames(Additive_2gen) <- A_2gen$pedigree$PEDIGREE_NAME
  rownames(Additive_2gen) <- A_2gen$pedigree$PEDIGREE_NAME
  
  return(Additive_2gen)
}


A_2gen_fresh <- Amat_2gen(fresh_data)
hist(diag(A_2gen_fresh))
hist(A_2gen_fresh[lower.tri(A_2gen_fresh, diag = F)])

s3saveRDS(A_2gen_fresh, bucket = bucket, object = 'ycao1/ITG_data/SweetCorn/A_2gen_fresh_2.rds')

A_2gen_proc <- Amat_2gen(proc_data)
hist(diag(A_2gen_proc))
hist(A_2gen_proc[lower.tri(A_2gen_proc, diag = F)])

s3saveRDS(A_2gen_proc, bucket = bucket, object = 'ycao1/ITG_data/SweetCorn/A_2gen_proc_2.rds')

############################### Generate Deep Pedigree Matrix ##############################

source('Code/Shared/deepped.R')


Amat_5gen <- function(dat, n_gen) {
  
  if (!'P1' %in% colnames(dat)){
    print("parentage pedigree doesn't exit")
    dat <- dat %>% 
      separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
    p1 <- unlist(lapply(dat$P1, FUN = pedSingle))
    p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    dat$P1 <- p1
    dat$P2[dat$P2=='']=NA
    p2 <- unlist(lapply(dat$P2, FUN = pedSingle))
    p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    dat$P2 <- p2
  }
  raw_dat_sub <- dat %>% 
    mutate(P1 = as.character(P1), 
           P2 = as.character(P2),
           trait = OBSRVTN_REF_CD, 
           TRAIT_VALUE = as.numeric(TRAIT_VALUE), 
           GROWSEASON = as.character(GROWSEASON),
           pedigree = as.character(PEDIGREE_NAME))
  raw_dat_sub[is.na(raw_dat_sub$P1), 'P2'] <- NA
  raw_dat_sub[is.na(raw_dat_sub$P2), 'P1'] <- NA
  
  
  # ped_file <- raw_dat_sub %>% 
  #   dplyr::select(PEDIGREE_NAME, P1, P2) %>% 
  #   distinct_all()
  
  new_token(fingerprint = F)
  gen5 <- make_ped_file(raw_dat_sub, crop_ids = CropIDs, n_gen = n_gen)
  A_matrix5 <- asreml.Ainverse(gen5$ped, fgen = list("inbreeding", n_gen))
  
  A_inv_5gen <- A_matrix5$ginv
  A_inv_sparse_5gen <- asreml.sparse2mat(A_inv_5gen)
  Additive_5gen <- round(solve(A_inv_sparse_5gen),2) 
  
  # A_matrix5_1 <- left_join(A_matrix5$pedigree, gen5$distinct_lines, by  = 'ID')
  
  A_matrix5_1 <- inner_join(A_matrix5$pedigree, gen5$distinct_lines, by  = 'ID')
  colnames(Additive_5gen) <- A_matrix5$pedigree$ID
  rownames(Additive_5gen) <- A_matrix5$pedigree$ID
  
  Additive_5gen_sub <- Additive_5gen[rownames(Additive_5gen) %in% A_matrix5_1$ID, colnames(Additive_5gen) %in% A_matrix5_1$ID]
  A_matrix5_1 <- A_matrix5_1[match(rownames(Additive_5gen_sub), A_matrix5_1$ID), ]
  
  rownames(Additive_5gen_sub) <- A_matrix5_1$pedigree
  colnames(Additive_5gen_sub) <- A_matrix5_1$pedigree
  
  
  return(Additive_5gen_sub)
}

A_5gen_fresh <- Amat_5gen(fresh_data, 5)
s3saveRDS(Additive_5gen_sub, bucket = bucket, object = 'ycao1/ITG_data/SweetCorn/A_5gen_fresh.rds')
