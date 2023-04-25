# This Script to calculate  BLUP for sweet corn

library(readxl)
library(tidyverse)
library(aws.s3)
library(plyr)
library(dplyr)

# how to remove the outliers automatically?
# Row number can be used as unique identifier.

######################################## Clean Data ######################################
options(scipen = 999)
sheet_names <- excel_sheets('Corn/sweet_corn_outliers_final.xlsx')
source('Code/Shared/misFunctions.R')

AWS.Authorization('ycao1')

path <- 'ycao1/ITG_data/SweetCorn/Data/'
bucket <- 'genome-analytics-perm-space'

files <- get_bucket(bucket = bucket, prefix = paste(path, 'Raw/', sep = ""))

for (i in 1:length(sheet_names)){
  print(sheet_names[i])
  file_name_1 <- gsub('outliers', 'Raw', sheet_names[i])
  path_rawfile_1 <- paste(path, 'Raw/', file_name_1, sep = '')
  raw_dat <- read.csv(text = rawToChar(get_object(path_rawfile_1, bucket = bucket)))
  
  # file_name_2 <- gsub('outliers', 'Raw_ELEN', sheet_names[i])
  # path_rawfile_2 <- paste(path, 'Raw/', file_name_2, sep = '')
  # raw_dat_2 <- read.csv(text = rawToChar(get_object(path_rawfile_2, bucket = bucket)))
  
  # raw_dat <- rbind(raw_dat_1, raw_dat_2)
  
  outlier_file <- read_excel('Corn/sweet_corn_outliers_final.xlsx', sheet = sheet_names[i])
  outlier_row <- outlier_file$Row[grep('Discard', outlier_file$`Decision - Breeder`, ignore.case = TRUE)]
  
  cln_dat <- raw_dat %>% 
    filter(!V1 %in% outlier_row)
  s3write_using(cln_dat,
                FUN = write.csv,
                bucket = bucket,
                object = paste(path, 'Raw/', gsub('outliers', 'clean', sheet_names[i]), sep = ''))
}


########################################## BLUP Calculation #############################################
library(asreml)
library(asremlPlus)
library(doBy)

Crop <- 'Corn'
source('Code/Shared/deepped.R')
source('Code/Shared/BLUPFunctions.R')
# PCM0_sub <- PCM0_dat %>% 
#   filter(!GROWSEASON %in% c('2019:04'))

metric <- 'ABLUP'

fresh_data <- data.frame()
proc_data <- data.frame()

file_list <- get_bucket(bucket = bucket, prefix = paste(path, 'Raw/', sep = ""))
file_list <- file_list[grepl('_Raw', file_list)]

AWS.Authorization('ycao1')
for (i in 1:length(file_list)){
  print(file_list[[i]]$Key)
  infile <- paste0('s3://','genome-analytics-perm-space/', file_list[[i]]$Key)
  csvfile <- get_object(infile)
  csvfile_1 <- rawToChar(csvfile)
  con <- textConnection(csvfile_1)
  temp_dat <- read.csv(con)
  close(con)
  if (grepl('_Fresh_', file_list[[i]]$Key)){
    fresh_data<- plyr::rbind.fill(fresh_data, temp_dat)
  }else{
    proc_data <- plyr::rbind.fill(proc_data, temp_dat)
  }
}

start_time <- Sys.time()
if (metric == 'ABLUP'){
  fresh_blup <- Ablup_fun(fresh_data, n_gen = 5, CropIDs = CropIDs)
  psca_fresh <- fresh_blup[[2]] %>% 
    filter(!is.na(PEDIGREE_NAME))
  pgca_fresh <- fresh_blup[[1]] %>% 
    filter(!is.na(PEDIGREE_NAME))
  
  s3write_using(psca_fresh, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/PSCA_Fresh_2014_2018_deep.csv')
  
  s3write_using(pgca_fresh, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/PGCA_Fresh_2014_2018_deep.csv')
  
  proc_blup <- Ablup_fun(proc_data, n_gen = 5, CropIDs = CropIDs)
  psca_proc <- proc_blup[[2]] %>% 
    filter(!is.na(PEDIGREE_NAME))
  pgca_proc <- proc_blup[[1]] %>% 
    filter(!is.na(PEDIGREE_NAME))
  
  s3write_using(psca_proc, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/PSCA_Proc_2014_2018_deep.csv')
  
  s3write_using(pgca_proc, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/PGCA_Proc_2014_2018_deep.csv')
  
  # summary data: heritability, N, n_hybrids
  
  fresh_h2 <- unique(psca_fresh[,c('trait', 'h2', 'N')])
  n_hybrids <- as.data.frame(table(psca_fresh$trait))
  colnames(n_hybrids) <- c('trait', 'unique_hybrids')
  fresh_h2 <- fresh_h2 %>% 
    join(n_hybrids, by = 'trait')
  s3write_using(fresh_h2, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/Fresh_heritability_deep.csv')
  
  proc_h2 <- unique(psca_proc[,c('trait', 'h2', 'N')])
  n_hybrids <- as.data.frame(table(psca_proc$trait))
  colnames(n_hybrids) <- c('trait', 'unique_hybrids')
  proc_h2 <- proc_h2 %>% 
    join(n_hybrids, by = 'trait')
  s3write_using(proc_h2, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/Proc_heritability_deep.csv')
  
}else{
  fresh_blup <- pblup_fun(fresh_data)
  psca_fresh <- fresh_blup[[2]]
  pgca_fresh <- fresh_blup[[1]]
  
  s3write_using(psca_fresh, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/PSCA_Fresh_2014_2018.csv')
  
  s3write_using(pgca_fresh, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/PGCA_Fresh_2014_2018.csv')
  
  proc_blup <- pblup_fun(proc_data)
  psca_proc <- proc_blup[[2]]
  pgca_proc <- proc_blup[[1]]
  
  s3write_using(psca_proc, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/PSCA_Proc_2014_2018.csv')
  
  s3write_using(pgca_proc, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/PGCA_Proc_2014_2018.csv')
  
  # summary data: heritability, N, n_hybrids
  
  fresh_h2 <- unique(psca_fresh[,c('trait', 'h2', 'N')])
  n_hybrids <- as.data.frame(table(psca_fresh$trait))
  colnames(n_hybrids) <- c('trait', 'unique_hybrids')
  fresh_h2 <- fresh_h2 %>% 
    join(n_hybrids, by = 'trait')
  s3write_using(fresh_h2, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/Fresh_heritability.csv')
  
  proc_h2 <- unique(psca_proc[,c('trait', 'h2', 'N')])
  n_hybrids <- as.data.frame(table(psca_proc$trait))
  colnames(n_hybrids) <- c('trait', 'unique_hybrids')
  proc_h2 <- proc_h2 %>% 
    join(n_hybrids, by = 'trait')
  s3write_using(proc_h2, FUN = write.csv, 
                bucket = bucket,
                object = 'ycao1/ITG_data/SweetCorn/BLUP/Proc_heritability.csv')
}

Sys.time() - start_time
