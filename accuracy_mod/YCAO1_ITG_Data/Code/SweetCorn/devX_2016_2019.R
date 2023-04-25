# Find dev crosses within a year 


library(aws.s3)
library(tidyverse)
library(caret)
library(glmnet)
library(jsonlite)
library(data.table)
library(Rcpp)
library(httr)


source('vaultCredentials.R')
# source('Code/Shared/ElasticNet.R')

crop <- 'Corn'
bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'


VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")


# Get Training data

############################ Save Training Data in S3 ############################

s3load(object = 'Corn/Corn_Parentals.RData', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')

years <- c('2016', '2017', '2018','2019', '2020')
AllData <- data.frame()
for (selectedYear in years){
  infile <- paste0(path,crop,"/",selectedYear,"/AllYearData.RData")
  s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")
  AllData <- plyr::rbind.fill(AllData,YearPhenos)
}

seasons <- c('2016:10', '2017:04', '2017:10', '2018:04', '2018:10', '2019:04', '2019:10', '2020:04')

AllData_sub <- AllData %>% 
  filter(LTYPE == 'SegPop', 
         TRIAL_TYPE %in% c('Nursery', 'YieldTrial')) %>% 
  mutate(devX_2019 = grepl('19 10|20 04', substring(SOURCE, 1, 5)),
         devX_2018 = grepl('18 10|19 04', substring(SOURCE, 1, 5)),
         devX_2017 = grepl('17 10|18 04', substring(SOURCE, 1, 5)),
         devX_2016 = grepl('16 10|17 04', substring(SOURCE, 1, 5))) %>% 
  mutate(devX = devX_2019 + devX_2018 + devX_2017 + devX_2016) %>% 
  filter(devX != 0)


AllData_sub_1 <- AllData_sub %>% 
  inner_join(CropParents, by =  c('PEDIGREE_NAME' = 'ChildPedigree'))


AllData_sub_2 <- AllData_sub_1 %>% 
  select(PEDIGREE_NAME, SOURCE, P1, P2, LTYPE, ChildGeneration, GROWSEASON, devX_2016, devX_2017, devX_2018, devX_2019, devX) %>% 
  filter(ChildGeneration == 'F1') %>% 
  unique() %>% 
  rename(GENERATION = ChildGeneration)


write.csv(AllData_sub_2, 'ProcTomato/DevX_2016_2019.csv', row.names = F)

