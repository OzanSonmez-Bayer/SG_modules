# BLUP calculation
library(asreml)
library(asremlPlus)
library(aws.s3)
library(plyr)
library(dplyr)
library(doBy)
library(readxl)
library(stringr)
options(scipen = 999)

source('Code/Shared/misFunctions.R')
path <- ('ycao1/ITG_data/ProcTomato/')
bucketname <- 'genome-analytics-perm-space'

AWS.Authorization('ycao1')
# source('Credential/Credential_yc.R')

Crop <- 'Tomato'
# Get a list of files in the path
dfBucket <- get_bucket_df(bucketname, paste0(path, 'Data'))
file_list <- dfBucket$Key[grep('PCM', dfBucket$Key)]
file_PCM0 <- file_list[grep('PCM0', file_list)]
file_PCM1 <- file_list[grep('PCM1', file_list)]

# Combined PCM 0 data
PCM0_dat <- data.frame()
for (i in 1:length(file_PCM0)){
  print(file_PCM0[i])
  infile <- paste0('s3://','genome-analytics-perm-space/', file_PCM0[i])
  csvfile <- get_object(infile)
  csvfile_1 <- rawToChar(csvfile)
  con <- textConnection(csvfile_1)
  temp <- read.csv(con)
  close(con)
  if(grep(file_PCM0[i], pattern = 'PCM0')>=1){
    if (i == 1){
      PCM0_dat <- temp
    }else{
      PCM0_dat <- rbind.fill(PCM0_dat, temp)
    }
  }
}
# Combine PCM 1 data
PCM1_dat <- data.frame()
for (j in 1:length(file_PCM1)){
  print(file_PCM1[j])
  infile <- paste0('s3://','genome-analytics-perm-space/', file_PCM1[j])
  csvfile <- get_object(infile)
  csvfile_1 <- rawToChar(csvfile)
  con <- textConnection(csvfile_1)
  temp <- read.csv(con)
  close(con)
  if (grep(file_PCM1[j], pattern = 'PCM1')>=1){
    if (j == 1){
      PCM1_dat <- temp
    }else{
      PCM1_dat <- rbind.fill(PCM1_dat, temp)
    }
  }
}

# Check UOM 
PCM0_dat <- conv_uom(PCM0_dat)
PCM1_dat <- conv_uom(PCM1_dat)
# Test code
source('Code/Shared/deepped.R')
source('Code/Shared/BLUPFunctions.R')
# PCM0_sub <- PCM0_dat %>% 
#   filter(!GROWSEASON %in% c('2019:04'))

metric <- 'ABLUP'

if (metric == 'ABLUP'){
  # PBLUP with A matrix built based on deeper pedigree. 
  A_blup <- Ablup_fun(PCM0_dat, n_gen = 5, CropIDs = CropIDs)
  
  pgca_pcm0 <- A_blup[[1]]
  psca_pcm0 <- A_blup[[2]]
  s3write_using(pgca_pcm0, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'PGCA_PCM0_2014_2019_deep.csv'))
  s3write_using(psca_pcm0, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'PSCA_PCM0_2014_2019_deep.csv'))
  
  
  # PCM1_sub <- PCM1_dat %>% 
  #   filter(!GROWSEASON %in% c('2019:04'))
  
  blup_pcm1 <- Ablup_fun(PCM1_dat, n_gen = 5, CropIDs = CropIDs)
  pgca_pcm1 <- blup_pcm1[[1]]
  psca_pcm1 <- blup_pcm1[[2]]
  s3write_using(pgca_pcm1, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'PGCA_PCM1_2014_2019_deep.csv'))
  s3write_using(psca_pcm1, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'PSCA_PCM1_2014_2019_deep.csv'))
  # Trait Summary
  # dat <- read.csv('ProcTomato/PSCA_PCM0_2014_2018.csv')
  # h2_psca_pcm0 <- unique(dat[,c('trait', 'h2', 'N')])
  # h2_psca_pcm0$n_hybrids <- unname(table(dat$trait)[h2_psca_pcm0$trait])
  
  pcm0_h2 <- unique(psca_pcm0[,c('trait', 'h2', 'N')])
  n_hybrids <- as.data.frame(table(psca_pcm0$trait))
  colnames(n_hybrids) <- c('trait', 'unique_hybrids')
  pcm0_h2 <- pcm0_h2 %>% 
    join(n_hybrids, by = 'trait')
  s3write_using(pcm0_h2, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'h2_PCM0_deep.csv'))
  
  
  
  pcm1_h2 <- unique(psca_pcm1[,c('trait', 'h2', 'N')])
  n_hybrids <- as.data.frame(table(psca_pcm1$trait))
  colnames(n_hybrids) <- c('trait', 'unique_hybrids')
  pcm1_h2 <- pcm1_h2 %>% 
    join(n_hybrids, by = 'trait')
  s3write_using(pcm1_h2, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'h2_PCM1_deep.csv'))
  
}else{
  
  tst_blup <- pblup_fun(PCM0_dat)
  
  pgca_pcm0 <- tst_blup[[1]]
  psca_pcm0 <- tst_blup[[2]]
  s3write_using(pgca_pcm0, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'PGCA_PCM0_2014_2019.csv'))
  s3write_using(psca_pcm0, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'PSCA_PCM0_2014_2019.csv'))
  
  
  # PCM1_sub <- PCM1_dat %>% 
  #   filter(!GROWSEASON %in% c('2019:04'))
  
  blup_pcm1 <- pblup_fun(PCM1_dat)
  pgca_pcm1 <- blup_pcm1[[1]]
  psca_pcm1 <- blup_pcm1[[2]]
  s3write_using(pgca_pcm1, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'PGCA_PCM1_2014_2019.csv'))
  s3write_using(psca_pcm1, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'PSCA_PCM1_2014_2019.csv'))
  # Trait Summary
  # dat <- read.csv('ProcTomato/PSCA_PCM0_2014_2018.csv')
  # h2_psca_pcm0 <- unique(dat[,c('trait', 'h2', 'N')])
  # h2_psca_pcm0$n_hybrids <- unname(table(dat$trait)[h2_psca_pcm0$trait])
  
  pcm0_h2 <- unique(psca_pcm0[,c('trait', 'h2', 'N')])
  n_hybrids <- as.data.frame(table(psca_pcm0$trait))
  colnames(n_hybrids) <- c('trait', 'unique_hybrids')
  pcm0_h2 <- pcm0_h2 %>% 
    join(n_hybrids, by = 'trait')
  s3write_using(pcm0_h2, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'h2_PCM0.csv'))
  
  
  
  pcm1_h2 <- unique(psca_pcm1[,c('trait', 'h2', 'N')])
  n_hybrids <- as.data.frame(table(psca_pcm1$trait))
  colnames(n_hybrids) <- c('trait', 'unique_hybrids')
  pcm1_h2 <- pcm1_h2 %>% 
    join(n_hybrids, by = 'trait')
  s3write_using(pcm1_h2, FUN = write.csv, 
                bucket = bucketname,
                object = paste0(path, 'BLUP/', 'h2_PCM1.csv'))
}


################ Generate A matrices ############
source('Code/Shared/AmatCalc.R')

A_5gen <- Amat_5gen(PCM0_dat, n_gen = 5)

A_2gen <- Amat_2gen(PCM0_dat)
