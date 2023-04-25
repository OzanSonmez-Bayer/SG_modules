# Pull Sweet Corn Data

library(aws.s3)
library(plyr)
library(tidyverse)
library(readxl)

############################### Check Entry ##############################
set_df <- read_excel('Corn/Screening sets 2014-2018.xlsx')

source('Credential/Credential_h2h.R')

path <- '/shared/Veg_Phenotypes/H2H_Data/'
crop <- 'Corn'
selectedProg <- '6S'

raw_data <- data.frame()
years <- unique(set_df$Season)

for (SelectedSeason in years){
  sets <- set_df$TEST_SET_NAME[which(set_df$Season == SelectedSeason)]
  SelectedSeason <- gsub(':', '/', SelectedSeason)
  infile <- paste0(path,crop,"/",SelectedSeason,"/",selectedProg,".RData")
  
  if(head_object(object = infile, bucket = "genome-analytics-perm-space")[[1]]){
    print("OBJECT FOUND")
    s3load(object = infile, bucket = "genome-analytics-perm-space")
    Phenos <- Phenos %>% 
      filter(TEST_SET_NAME %in% sets)
    raw_data <- plyr::rbind.fill(raw_data, Phenos)
  } else {
    print("OBJECT NOT FOUND")
  }
}

# Check Entry number (Unique pedigrees) by TEST_SET_NAME
set_df <- set_df %>% 
  select(Season, TEST_SET_NAME, Type, `# Entries (P/T)`) %>% 
  `colnames<-`(c('GROWSEASON', 'TEST_SET_NAME', 'Type', 'EntryNumber'))

raw_data <- raw_data %>% 
  mutate(Type = NA) %>% 
  mutate(Type = if_else(TEST_SET_NAME %in% set_df$TEST_SET_NAME[which(set_df$Type == 'Fresh')], 'Fresh', 'Processing'))

entry_num <- aggregate(data = raw_data, PEDIGREE_NAME ~ GROWSEASON + TEST_SET_NAME + Type, function(x) length(unique(x)))

# Compare the entry number with the one from breeder

compare_entry <- merge(entry_num, set_df, by = c('GROWSEASON', 'TEST_SET_NAME', 'Type'))

write.csv(compare_entry, 'Corn/check_entry.csv') # We got all most all the entries.


################################## Outlier Analysis ############################
# Load data with traits of interest

path_corn<- 'ycao1/ITG_data/SweetCorn/Data/'
bucket_name <- 'genome-analytics-perm-space'
source('Credential_yc.R')

install.packages('xlsx')
library(xlsx)

files <- get_bucket(bucket = bucket_name, prefix = paste(path_corn, 'Outliers/', sep = ""))

tempfile <- read.csv(text = rawToChar(get_object(files[[2]]$Key, bucket = bucket_name)))
xlsx::write.xlsx(tempfile,file='Corn/sweet_corn_outliers.xlsx',sheetName = strsplit(files[[2]]$Key, split = '/')[[1]][6])

for (i in 3:length(files)){
  print(files[[i]]$Key)
  tempfile <- read.csv(text = rawToChar(get_object(files[[i]]$Key, bucket = bucket_name)))
  xlsx::write.xlsx(tempfile,file='Corn/sweet_corn_outliers.xlsx',sheetName = strsplit(files[[i]]$Key, split = '/')[[1]][6], append = TRUE)
}


# Get mean and standard deviation for each observatons
outDF_fresh <- data.frame()
outDF_proc <- data.frame()

for(i in 2:length(files)){
  print(files[[i]]$Key)
  tempfile <- read.csv(text = rawToChar(get_object(files[[i]]$Key, bucket = bucket_name)))
  file_name <- gsub('outliers', 'Raw', strsplit(files[[i]]$Key, split = '/')[[1]][6])
  rawfile_name <- paste(path_corn, 'Raw/', file_name, sep = '')
  print(rawfile_name)
  temp_raw <- read.csv(text = rawToChar(get_object(rawfile_name, bucket = bucket_name)))
  
  tempfile$condition <- paste(tempfile$Cross, tempfile$SetName, tempfile$Season, tempfile$Trait, sep = "__")
  dat <- temp_raw %>% 
    select(V1, CROSS_NAME, TEST_SET_NAME, PLOT_BID,GROWSEASON, OBSRVTN_REF_CD, TRAIT_VALUE) %>%
    mutate(condition = paste(CROSS_NAME, TEST_SET_NAME, GROWSEASON, OBSRVTN_REF_CD, sep = "__")) %>% 
    filter(condition %in% tempfile$condition)
  
  dat_summary <- dat  %>% 
    group_by(condition) %>% 
    summarise(avg = round(mean(as.numeric(TRAIT_VALUE, na.rm = T)),2), 
              median = round(median(as.numeric(TRAIT_VALUE, na.rm = T)),2),
              sd = round(sd(as.numeric(TRAIT_VALUE, na.rm = T)),2)) %>% 
    separate(condition, c('CROSS_NAME', 'TEST_SET_NAME', 'GROWSEASON', 'OBSRVTN_REF_CD'),sep = "__", remove = TRUE) %>% 
    right_join(dat, by = c('CROSS_NAME', 'TEST_SET_NAME', 'GROWSEASON', 'OBSRVTN_REF_CD')) %>% 
    select(-condition)
  dat_summary$outlier <- ifelse(dat_summary$V1 %in% tempfile$Row, 1, 0) # 1 is a detected outlier. 
  
  if (grepl('Fresh', files[[i]]$Key)){
    if(i == 2){
      outDF_fresh <- dat_summary
    }else{
      outDF_fresh <- rbind(outDF_fresh, dat_summary)
    }
  }else{
    if (i == 2){
      outDF_proc <- dat_summary
    }else{
      outDF_proc <- rbind(outDF_proc, dat_summary) 
    }
  }
}

write.csv(outDF_proc, 'Proc_outliers_obs.csv', row.names = F)
write.csv(outDF_fresh, 'Fresh_outliers_obs.csv', row.names = F)