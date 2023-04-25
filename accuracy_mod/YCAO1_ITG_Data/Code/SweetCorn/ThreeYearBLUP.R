# Three-year BLUP

library(asreml)
library(aws.s3)
library(tidyverse)
library(dplyr)
library(httr)
library(asremlPlus)
# install.packages('/mnt/packages/doBy_4.6-2.tar.gz',repos = NULL, type = 'source')
library(doBy)
library(tidyr)


bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
Crop <- 'Corn'

# Credentials

source('vaultCredentials.R')
source('Code/Shared/BLUPFunctions.R')
source('Code/Shared/deepped.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

# Load Crop IDs
# s3load(object = 'Corn/Corn_IDs.RData', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')

# Pull phenotypic data from H2H_Data
years <- c('2018','2019','2020')
AllData <- data.frame()
for (selectedYear in years){
  infile <- paste0(path,Crop,"/",selectedYear,"/AllYearData.RData")
  s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")
  AllData <- rbind(AllData,YearPhenos)
}


seasons <- c('2018:04', '2019:04', '2020:04')
pattern <- c('1E', '2E', '3E', '1F', '2F', '3F')
stage <- c('S1', 'S2', 'S3')
# selectedSets <- grep(paste('^',paste(pattern,collapse="|^"), sep = ''), unique(AllData$TEST_SET_NAME), value = TRUE)

traits <- c('ASI50D', 'EDIA', 'EHT','ELEN', 'ERSC', 'ESCCA', 'ESHK', 'FINAL','HSC','P50_BE','PA','PHT',
            'QUAL', 'RTLR','S50_BE', 'SC_ELEN', 'SC_FR', 'SC_HE', 'SC_RWCT', 'SC_SE', 'SC_TF', 'SCROW', 'SDV','TILLR')

# fresh_sc123 <- AllData %>% 
#   filter(GROWSEASON %in% seasons & EXPER_STAGE_REF_ID %in% c('S1', 'S2', 'S3') & 
#            TEST_SET_NAME %in% selectedSets & TRIAL_TYPE == 'YieldTrial' &
#            OBSRVTN_REF_CD %in% traits) %>% 
#   mutate(TREP = paste(REP_NUMBER, TRACK_ID, sep = '_'))

getpheno <- function(dat, sets_pattern, stage){
  selectedSets <- grep(paste('^',paste(sets_pattern,collapse="|^"), sep = ''), unique(AllData$TEST_SET_NAME), value = TRUE)
  bdat <- dat %>% 
    filter(GROWSEASON %in% seasons & EXPER_STAGE_REF_ID %in% stage  & 
             TEST_SET_NAME %in% selectedSets & TRIAL_TYPE == 'YieldTrial' &
             OBSRVTN_REF_CD %in% traits) %>% 
    mutate(TREP = paste(REP_NUMBER, TRACK_ID, sep = '_'))
}

long_to_wide <- function(dat){
  dat_wide <- data.frame()
  for (trt in 1: length(unique(dat$trait))){
    temp <- dat %>% 
      dplyr::filter(trait == unique(dat$trait)[trt])
    colnames(temp)[-1] <- paste(unique(dat$trait)[trt],colnames(fresh_sc123_blup[[2]])[-1], sep = '_')
    
    if(nrow(dat_wide) == 0){
      dat_wide <- temp
    }else{
      dat_wide <- dat_wide %>% 
        full_join(temp, by = 'PEDIGREE_NAME')
    }
  }
  return(dat_wide)
}


fresh_sc123 <- getpheno(AllData, sets_pattern = pattern, stage = stage)

# BLUP
fresh_sc123_blup <- Ablup_fun(dataset = fresh_sc123, n_gen = 5, CropIDs = CropIDs) 

# Convert to wide format

fresh_sc123_sca <- long_to_wide(fresh_sc123_blup[[2]])
fresh_sc123_gca <- long_to_wide(fresh_sc123_blup[[1]])

write.csv(fresh_sc123_sca, file = 'Corn/fresh_sc123_sca.csv', row.names = F)
write.csv(fresh_sc123_gca, file = 'Corn/fresh_sc123_gca.csv', row.names = F)


#### Add GCV HETGP 

keyline <- readxl::read_excel('/domino/datasets/local/SweetCorn/Key Inbred File for Yujing.xlsx', sheet = '2019 key inbred P1+')
keyline <- keyline %>% 
  select(-Pedigree...27) %>% 
  rename(Pedigree = Pedigree...62)

keyline <- keyline[,1:61]

blup <- s3read_using(FUN = read.csv, 
                     object = 'ycao1/OriginSelection/Corn/ThreeYear_BLUP/GCA/Temperate_ThreeYear_2019.csv', 
                     bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

blup_keyline <- keyline %>% 
  left_join(blup, by = c('Pedigree' = 'PEDIGREE_NAME'))

write.csv(blup_keyline, file = '/mnt/Temperate_ThreeYear_2019_keyline.csv', row.names = F)

# s3write_using(blup_keyline, FUN = write_csv, 
#               object = 'ycao1/OriginSelection/Corn/ThreeYear_BLUP/GCA/Temperate_ThreeYear_2016_keyline.csv',
#               bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
