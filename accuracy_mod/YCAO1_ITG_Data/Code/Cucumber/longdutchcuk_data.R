# Pull Sweet Corn Data

library(aws.s3)
library(readxl)
library(plyr)
library(dplyr)
library(tidyverse)
library(outliers)
library(xlsx)

############################# Check Entry ###############################
if (grepl('WINDOWS', Sys.getenv('windir'))){
  setwd('C:/Users/YCAO1/Desktop/P-1/VEG_ITG/Cucumber')
}

source('Code/Shared/misFunctions.R')
set_df <- read.csv('Cucumber/set_df.csv')
traits <- read_excel('Cucumber/Y3Traits.xlsx')

traits_pcm1 <- traits$FtsTrait[which(traits$PCM1 == 1)]
traits_pcm2 <- traits$FtsTrait[which(traits$PCM2 == 1)]


crop <- 'Cucumber'
selectedProg <- 'Y3'

raw_data <- data.frame()
years <- unique(set_df$season)[!is.na(unique(set_df$season))]
for (SelectedSeason in years){
  sets <- unique(set_df$set.name[which(set_df$season == SelectedSeason)])
  temp_data <- pullData(season = SelectedSeason, 
                        setlist = sets, 
                        trtlist = traits_pcm1, 
                        crop = crop, 
                        selectedProg = selectedProg)
  
  raw_data <- rbind.fill(raw_data, temp_data)
}

# Check UOM
raw_data <- conv_uom(raw_data)

# Check Entry Number by Test Set Name
# PCM 1
entry_num_pcm1 <- aggregate(data = raw_data[raw_data$OBSRVTN_REF_CD %in% traits_pcm1, ], 
                       PEDIGREE_NAME ~ GROWSEASON + TEST_SET_NAME, function(x) length(unique(x)))

# Check Breeder Provided Info

breeder <- data.frame()
for (i in c(2:6)){
  temp_dat <- read_excel('Cucumber/Y3_TestSets.xlsx', sheet = i)
  temp_agg <- aggregate(data = temp_dat, `Pedigree HYBR` ~ season + `set name`, function(x)length(unique(x)))
  breeder <- rbind(breeder, temp_agg)
}

# Compare raw data with breeder's info

compare_entry <- breeder %>% 
  `colnames<-`(c('GROWSEASON', 'TEST_SET_NAME', 'PEDIGREE_NAME')) %>% 
  full_join(entry_num_pcm1, by = c('GROWSEASON', 'TEST_SET_NAME'))

total_entry <- compare_entry %>% 
  group_by(GROWSEASON) %>% 
  summarise(total_entry = sum(PEDIGREE_NAME.y, na.rm = T))

# Data from 2018:05 is missing from raw data: 15suv29, 15suv31, 15SUV41, and 15CMK

# GROWSEASON total_entry
# 2015:01            488
# 2015:05            412
# 2015:08            496
# 2016:01            478
# 2016:05            488
# 2016:08            472
# 2017:01            468
# 2017:05            467
# 2017:08            456
# 2018:01            466
# 2018:05            415
# 2018:08            416
# 2019:01            556
# 2019:05            602
# 2019:08            576

# Check number of observations

total_obs <- raw_data %>% 
  filter(OBSRVTN_REF_CD %in% traits_pcm1) %>% 
  group_by(GROWSEASON) %>% 
  count()
  
# GROWSEASON    n
# 2015:01    141764
# 2015:05    191840
# 2015:08    198961
# 2016:01    163867
# 2016:05    286302
# 2016:08    187885
# 2017:01    173179
# 2017:05    327829
# 2017:08    252990
# 2018:01    272187
# 2018:05    246555
# 2018:08    155805
# 2019:01    249232
# 2019:05    235202
# 2019:08    312338

### Create a new for Cucumber in S3 to store results
source('Credential/Credential_yc.R')
bucket <- 'genome-analytics-perm-space'
path <- '/ycao1/ITG_data/Cucumber'
put_folder('/ycao1/ITG_data/Cucumber', bucket = bucket)
put_folder('/ycao1/ITG_data/Cucumber/Raw', bucket = bucket)

############################ Calculated Traits #############################

# Calculated traits are calculated wihtin each grow season
pcm1_dat <- raw_data
pcm1_dat_2 <- data.frame()
for (season in unique(pcm1_dat$GROWSEASON)){
  print(season)
  calc_dat <- pcm1_dat %>% 
    filter(GROWSEASON == season)
  
  calc_dat_1 <-cal_trt(calc_dat)
  pcm1_dat_2 <- rbind(pcm1_dat_2, calc_dat_1)
}

s3save(pcm1_dat_2, 
       object = '/ycao1/ITG_data/Cucumber/Raw/pcm1_dat.Rdata',
       bucket = bucket)

# Note: how to differentiate data from PCM 1 and PCM 2?
pcm2_dat <- pcm1_dat_2 %>% 
  filter(OBSRVTN_REF_CD %in% traits_pcm2)
# pcm2_dat_2 <- data.frame()
# for (season in unique(pcm2_dat$GROWSEASON)){
#   print(season)
#   calc_dat <- pcm2_dat %>% 
#     filter(GROWSEASON == season)
#   
#   calc_dat_1 <- cal_trt(calc_dat)
#   pcm2_dat_2 <- rbind(pcm2_dat_2, calc_dat_1)
# }
# s3write_using(pcm2_dat_2, 
#               FUN = write.csv, 
#               bucket = bucket,
#               object = paste(path, '/Raw/pcm2_dat.csv'), sep = '')
s3save(pcm2_dat, 
       object = '/ycao1/ITG_data/Cucumber/Raw/pcm2dat.Rdata',
       bucket = bucket)
################################ Outlier Analysis ######################

for (season in unique(pcm1_dat_2$GROWSEASON)){
  print(season)
  temp_out <- rm.out(pcm1_dat_2[which(pcm1_dat_2$GROWSEASON == season), ])
  
  xlsx::write.xlsx(temp_out, file = 'Cucumber/cucumber_outliers_pcm1.xlsx', 
                   sheetName = paste(unlist(strsplit(season, ":")), collapse = "_"), 
                   append = TRUE)
}

for (season in unique(pcm2_dat_2$GROWSEASON)){
  print(season)
  temp_out <- rm.out(pcm2_dat_2[which(pcm2_dat_2$GROWSEASON == season), ])
  
  xlsx::write.xlsx(temp_out, file = 'Cucumber/cucumber_outliers_pcm2.xlsx', 
                   sheetName = paste(unlist(strsplit(season, ":")), collapse = "_"), 
                   append = TRUE)
}
