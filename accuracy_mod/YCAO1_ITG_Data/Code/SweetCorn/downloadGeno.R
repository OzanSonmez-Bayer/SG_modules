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


# data <- s3read_using(read.table, object = "s3://genome-analytics-perm-space//shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/ProcessingTomato20200615-1_HDgenotypes.txt.gz")
# colnames(data) <- data[1,]

geno <- s3read_using(fread, 
                     object = "Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/HDgenotypes.txt.gz", 
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

meta <- s3read_using(fread, 
                     object = "Veg_Genotypes/Sweetcorn/firstpilot/SweetCorn_6S_ProgenyTable.txt",
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

markerdata <- s3read_using(fread, 
                           object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/MarkerMetadata.txt',
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

# source('Credential/Credential_sdi.R')
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
  dplyr::inner_join(geno_wide_012, by = 'ProgenyGermID')

s3write_using(geno_wide_ped, FUN = writing_csv, 
              object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/Geno_wide_SC1.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

# geno_wide_ped_1 <- geno_wide_ped[!duplicated(geno_wide_ped$ProgenyGermID),]


# AWS.Authorization('ycao1')
# s3write_using(geno_wide_ped, FUN = writing_csv, 
#               object = 'ycao1/ITG_data/SweetCorn/Geno/imputedGeno.csv', 
#               bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
# s3write_using(geno_wide_ped_1, FUN = write_csv, 
#               object = '/shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/imputedGeno_wide.csv', 
#               bucket = bucket)


samples_taqman <- as.numeric(unique(geno$ProgenyGermID[which(geno$ImputationStatus == 'fromLowDens')])) # 2720
samples_FP <- as.numeric(unique(geno$ProgenyGermID[which(geno$ImputationStatus == 'FPphased')])) # 582

ped_taqman <- geno_wide_ped[which(geno_wide_ped$ProgenyGermID %in% samples_taqman), c('ProgenyGermID', 'Pedigree')]
ped_FP <- geno_wide_ped[which(geno_wide_ped$ProgenyGermID %in% samples_FP), c('ProgenyGermID', 'Pedigree')]

s3write_using(ped_FP, FUN = writing_csv, 
              object = 'ycao1/ITG_data/SweetCorn/Geno/ped_FP.csv', 
              bucket = bucket)
s3write_using(ped_taqman, FUN = writing_csv, 
              object = 'ycao1/ITG_data/SweetCorn/Geno/ped_taqman.csv', 
              bucket = bucket)



####### Hybrids Geno ############

hybrid_geno <- s3read_using(fread, object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/HybridGenotypes/SweetCornImplementationHybridGT.txt.gz',
                            bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

hybrid_wide <- hybrid_geno %>% 
  dplyr::select(HybridGermID, MRN, Genotype) %>% 
  dcast(HybridGermID ~ MRN, value.var="Genotype")

# Turn MRN to 0, 1, 2
hybrid_wide_012 <- apply(hybrid_wide[, 2:ncol(hybrid_wide)], MARGIN = 2, 
                       FUN = function(x){
                         unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
                       })

hybrid_wide_012 <- cbind(hybrid_wide[,1], hybrid_wide_012)
hybrid_wide_ped <- hybrid_wide_012 %>% 
  mutate(HybridGermID = as.numeric(HybridGermID)) %>% 
  dplyr::inner_join(CropIDs_sub, by = c('HybridGermID' = 'ProgenyGermID'))

inbred_geno <- s3read_using(fread, object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/imputedGeno_wide.csv', 
                            bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

colnames(inbred_geno)[which(colnames(inbred_geno) == 'ProgenyGermID')] <- 'GermID'
colnames(hybrid_wide_ped)[which(colnames(hybrid_wide_ped) == 'HybridGermID')] <- 'GermID'


hybrid_wide_ped <- hybrid_wide_ped %>% 
  select(GermID, Pedigree, colnames(hybrid_wide_ped)[which(!colnames(hybrid_wide_ped) %in% c('GermID', 'Pedigree'))]) %>% 
  mutate(GermID = as.numeric(GermID))


s3write_using(hybrid_wide_ped, FUN = writing_csv, 
              object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/HybridGenotypes/HybridGeno_SC1_wide.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

s3write_using(geno_wide_ped, FUN = writing_csv, 
              object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/Geno_wide_SC1.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')


# Check to see if there are any NAs in hybrid gneo file


tst <- hybrid_wide_ped[, 3:ncol(hybrid_wide_ped)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
hybrid_wide_ped[,3:ncol(hybrid_wide_ped)] <- tst
dim(hybrid_wide_ped)

geno_1 <- s3read_using(fread, object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/imputedGeno_wide.csv', 
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

hybrid_geno_1 <- s3read_using(fread, object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/HybridGeno_wide.csv',
                              bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

tst <- hybrid_geno_1[, 3:ncol(hybrid_geno_1)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
hybrid_geno_1[,3:ncol(hybrid_geno_1)] <- tst

length(intersect(hybrid_geno_1$Pedigree, hybrid_wide_ped$Pedigree))
length(intersect(geno_1$Pedigree, geno_wide_ped$Pedigree))


hybrid_geno_1 <- hybrid_geno_1 %>% 
  filter(!Pedigree %in% intersect(hybrid_geno_1$Pedigree, hybrid_wide_ped$Pedigree))

geno_1 <- geno_1 %>% 
  filter(!Pedigree %in% intersect(geno_1$Pedigree, geno_wide_ped$Pedigree))


length(intersect(geno_wide_ped$Pedigree, hybrid_geno_1$Pedigree))

geno_combine <- do.call(rbind, list(geno_1, geno_wide_ped, hybrid_geno_1, hybrid_wide_ped))
str(geno_combine[,1:4])

s3saveRDS(geno_combine, object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds', 
          bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

