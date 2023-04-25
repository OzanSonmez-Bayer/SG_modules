library(aws.s3)
library(tidyverse)
library(caret)
library(glmnet)
library(jsonlite)
library(data.table)
library(Rcpp)
library(httr)


source('vaultCredentials.R')
source('Code/Shared/ElasticNet.R')

crop <- 'Tomato'
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
years <- c('2014', '2015','2016', '2017', '2018','2019')
AllData <- data.frame()
for (selectedYear in years){
  infile <- paste0(path,crop,"/",selectedYear,"/AllYearData.RData")
  s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")
  AllData <- plyr::rbind.fill(AllData,YearPhenos)
}

tst <- read.csv('ProcTomato/PCM0_data.csv')
seasons <- c('2014:04','2015:04', '2016:04', '2017:04', '2018:04', '2019:04')

setlists <- c(unique(as.character(tst$TEST_SET_NAME)), 'PCM1', 'PCM1TI', 'PMC1LATE', 'CKRT', 'AORT', 'AO_PCM2', 'CK_PCM2')
traits_itg <- c('FBRIX', 'AVGHBRIX', 'AVJB', 'FQUAL', 'FRFRMH', 'FZUNI', 'LBRIX', 'MAT', 'MATR', 'OST', 'PLCMP', 'PLTVG', 'PYLDPA', 'SAMPUW')

rm(tst)

itg_data <- AllData %>%
  filter(GROWSEASON %in% seasons & TEST_SET_NAME %in% setlists)

s3saveRDS(itg_data, object = 'ycao1/ITG/Tomato/9Z/pheno_training.rds',
          bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


########################### Format Geno ################################
geno <- s3read_using(fread, object = 'Tomato/ProcessingTomato/9Z/HybridGenotypes/ProcessingTomatoHybridGT.tab.gz', 
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

geno_wide <- geno %>% 
  rename(GermID = HybridGermID) %>% 
  select(GermID, MRN, Genotype) %>% 
  dcast(GermID ~ MRN, value.var="Genotype")

# Turn MRN to 0, 1, 2
geno_wide_012 <- apply(geno_wide[, 2:ncol(geno_wide)], MARGIN = 2, 
                       FUN = function(x){
                         unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
                       })

geno_wide_012 <- cbind(geno_wide[,1], geno_wide_012)

s3load(object = 'Tomato/Tomato_IDs.RData', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')

CropIDs_sub <- CropIDs %>% 
  dplyr::select(M.GERMPLASM.X_ID, M.GERMPLASM.PEDIGREE) %>% 
  dplyr::rename(GermID = M.GERMPLASM.X_ID, 
                Pedigree = M.GERMPLASM.PEDIGREE) %>% 
  mutate(GermID = as.numeric(GermID))

geno_wide_012$GermID <- as.numeric(geno_wide_012$GermID)

geno_wide_ped <- CropIDs_sub %>% 
  inner_join(geno_wide_012, by = 'GermID')

geno_wide_ped_1 <- geno_wide_ped[!duplicated(geno_wide_ped$GermID),]

# Fill NAs 

tst <- geno_wide_ped_1[, 3:ncol(geno_wide_ped_1)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
geno_wide_ped_1[,3:ncol(geno_wide_ped_1)] <- tst

s3saveRDS(geno_wide_ped_1, object = 'ycao1/ITG/Tomato/9Z/hybridGeno_training.rds',
          bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
