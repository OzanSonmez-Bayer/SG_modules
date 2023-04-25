# QC/QA

library(aws.s3)
library(tidyverse)

# Set working directory
setwd('C:/Users/ycao1/OneDrive - Monsanto/Migrated from My PC/Desktop/P-1/VEG_ITG/outliers')
prog <- '9Z'                 
crop <- 'Tomato'
path=('/shared/Veg_Phenotypes/H2H_Data/')


source('credential_h2h.R')

file_list_PCM0 <- list.files(path = 'PCM0' )
file_list_PCM1 <- list.files(path = 'PCM1')
outDF <- data.frame()

for (i in 1:length(file_list_PCM0)){
  print(file_list_PCM0[i])
  temp <- read.csv(paste0('PCM0/', file_list_PCM0[i], collapse = ''))
  season <- unique(temp$Season)
  season <- gsub(":","/",season)
  infile <- paste0(path,crop,"/",season,"/",prog,".RData")
  s3load(object = infile, bucket = "genome-analytics-perm-space")
  
  temp$condition <- paste(temp$Cross, temp$SetName, temp$Season, temp$Trait, sep = "/")
  dat <- Phenos %>% 
    select(CROSS_NAME, TEST_SET_NAME, GROWSEASON, OBSRVTN_REF_CD, TRAIT_VALUE) %>%
    mutate(condition = paste(CROSS_NAME, TEST_SET_NAME, GROWSEASON, OBSRVTN_REF_CD, sep = "/")) %>% 
    filter(condition %in% temp$condition)
  dat_summary <- dat  %>% 
    group_by(condition) %>% 
    summarise(avg = round(mean(as.numeric(TRAIT_VALUE)),2), 
              median = round(median(as.numeric(TRAIT_VALUE)),2),
              sd = round(sd(as.numeric(TRAIT_VALUE)),2)) %>% 
    separate(condition, c('CROSS_NAME', 'TEST_SET_NAME', 'GROWSEASON', 'OBSRVTN_REF_CD'),sep = "/", remove = TRUE) %>% 
    right_join(dat, by = c('CROSS_NAME', 'TEST_SET_NAME', 'GROWSEASON', 'OBSRVTN_REF_CD')) %>% 
    select(-condition)
  
  
  if (i == 1){
    outDF <- dat_summary
  }else{
    outDF <- rbind(outDF, dat_summary)
  }
}

write.csv(outDF, file = 'PCM0/all_obs_outlier.csv')

outDF <- data.frame()
for (i in 1:length(file_list_PCM1)){
  print(file_list_PCM1[i])
  temp <- read.csv(paste0('PCM1/', file_list_PCM1[i], collapse = ''))
  season <- unique(temp$Season)
  season <- gsub(":","/",season)
  infile <- paste0(path,crop,"/",season,"/",prog,".RData")
  s3load(object = infile, bucket = "genome-analytics-perm-space")
  
  temp$condition <- paste(temp$Cross, temp$SetName, temp$Season, temp$Trait, sep = "/")
  dat <- Phenos %>% 
    select(CROSS_NAME, TEST_SET_NAME, GROWSEASON, OBSRVTN_REF_CD, TRAIT_VALUE) %>%
    mutate(condition = paste(CROSS_NAME, TEST_SET_NAME, GROWSEASON, OBSRVTN_REF_CD, sep = "/")) %>% 
    filter(condition %in% temp$condition)
  dat_summary <- dat  %>% 
    group_by(condition) %>% 
    summarise(avg = round(mean(as.numeric(TRAIT_VALUE)),2), 
              median = round(median(as.numeric(TRAIT_VALUE)),2),
              sd = round(sd(as.numeric(TRAIT_VALUE)),2)) %>% 
    separate(condition, c('CROSS_NAME', 'TEST_SET_NAME', 'GROWSEASON', 'OBSRVTN_REF_CD'),sep = "/", remove = TRUE) %>% 
    right_join(dat, by = c('CROSS_NAME', 'TEST_SET_NAME', 'GROWSEASON', 'OBSRVTN_REF_CD')) %>% 
    select(-condition)
  
  
  if (i == 1){
    outDF <- dat_summary
  }else{
    outDF <- rbind(outDF, dat_summary)
  }
}


write.csv(outDF, file = 'PCM1/all_obs_outlier.csv')