# Check Processing Tomato Data

library(aws.s3)
library(plyr)
library(magrittr)
library(jsonlite)
library(data.table)
library(bit64)
library(tidyr)
library(httr)


source('Code/Shared/misFunctions.R')

options(scipen = 999)

bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-perm-space'
path <- '/H2H_Data/'

# Credentials

source('vaultCredentials.R')
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")


################################### Formatting Imputed Geno Files ############################

# Long format to wide format
# Add pedigree to the file 

infile <- 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/258linesFromThreeImportantOrigins_2021-07-27/HDgenotypes.txt.gz'
geno <- s3read_using(fread, 
                     object = infile, 
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

meta_infile <- 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/258linesFromThreeImportantOrigins_2021-07-27/SweetCorn_Temperate_300Lines_20210720_ProgenyTable.txt'
meta <- s3read_using(fread, 
                     object = meta_infile,
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

marker_infile <- "Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/258linesFromThreeImportantOrigins_2021-07-27/MarkerMetadata.txt"
markerdata <- s3read_using(fread, 
                           object = marker_infile,
                           bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

head(geno)
head(meta)
head(markerdata)

# Download reference to add pedigree to geno file
geno_wide <- geno %>% 
  dplyr::select(ProgenyGermID, MRN, HDgenotype) %>% 
  dcast(ProgenyGermID ~ MRN, value.var="HDgenotype")

# Turn MRN to 0, 1, 2
geno_wide_012 <- apply(geno_wide[, 2:ncol(geno_wide)], MARGIN = 2, 
                       FUN = function(x){
                         unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
                       })

geno_wide_012 <- cbind(geno_wide[,1], geno_wide_012)


# Pull Reference ID from S3

Crop <- 'Corn'

s3load(object = 'Corn/Corn_IDs.RData', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')

CropIDs_sub <- CropIDs %>%
  dplyr::select(M.GERMPLASM.X_ID,  M.GERMPLASM.PEDIGREE) %>%
  dplyr::rename(ProgenyGermID = M.GERMPLASM.X_ID,
                Pedigree = M.GERMPLASM.PEDIGREE) %>%
  mutate(ProgenyGermID = as.numeric(ProgenyGermID)) %>%
  unique()


geno_wide_012$ProgenyGermID <- as.numeric(geno_wide_012$ProgenyGermID)

geno_wide_ped <- CropIDs_sub %>% 
  dplyr::inner_join(geno_wide_012, by = 'ProgenyGermID') %>% 
  dplyr::rename(GermID = ProgenyGermID)

s3write_using(geno_wide_ped, FUN = writing_csv, 
              object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/258linesFromThreeImportantOrigins_2021-07-27/geno_wide_300lines.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

######################### Adding New Lines to Existing Imputed Geno File ################################

existing_geno <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds',
                           bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

existing_geno <- existing_geno %>% 
  dplyr::mutate(GermID = as.numeric(GermID))

keylines <- s3read_using(fread, 
                         object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/258linesFromThreeImportantOrigins_2021-07-27/germIDsOfThe258LinesForPrediction--REALLY.txt', 
                         bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

# Add Pedigree to the keylines 

keylines_ped <- keylines %>% 
  dplyr::rename(ProgenyGermID = V1) %>% 
  dplyr::mutate(ProgenyGermID = as.numeric(ProgenyGermID)) %>% 
  dplyr::left_join(CropIDs_sub, by = 'ProgenyGermID')

head(keylines_ped)

s3write_using(keylines_ped, FUN = writing_csv, 
              object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/258linesFromThreeImportantOrigins_2021-07-27/GermID_Ped_258lines_prediction.csv',
              bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

# Reference panel is included in the new line file. We need to check duplicated lines 

# It indicates overlapping lines, could be from reference panel
# We need clarification from GBD
# Note, a different reference panel is used for implementation. It is different from training set!

reference_panel <- intersect(geno_wide_ped$GermID, existing_geno$GermID) 

geno_wide_ped_sub <- geno_wide_ped %>% 
  dplyr::filter(!GermID %in% reference_panel)


combined_geno <- rbind(existing_geno, geno_wide_ped_sub)

s3saveRDS(combined_geno, object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all_07292021.rds', 
          bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')



