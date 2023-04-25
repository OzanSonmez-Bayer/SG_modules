# Check Processing Tomato Data

library(aws.s3)
library(plyr)
library(magrittr)
library(jsonlite)
library(data.table)
library(bit64)
library(tidyr)


source('Code/Shared/misFunctions.R')

options(scipen = 999)

AWS.Authorization('ycao1')

bucket <- 'genome-analytics-perm-space'
path <- '/shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/'


# data <- s3read_using(read.table, object = "s3://genome-analytics-perm-space//shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/ProcessingTomato20200615-1_HDgenotypes.txt.gz")
# colnames(data) <- data[1,]

geno <- s3read_using(fread, 
                     object = "s3://genome-analytics-perm-space//shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/ProcessingTomato20200615-1_HDgenotypes.txt.gz")

meta <- s3read_using(fread, 
                     object = "s3://genome-analytics-perm-space//shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/ProcessingTomato20200615-1_MarkerMetadata.txt")

head(geno)
head(meta)

# Download reference to add pedigree to geno file
geno_wide <- geno %>% 
  select(ProgenyGermID, MRN, HDgenotype) %>% 
  dcast(ProgenyGermID ~ MRN, value.var="HDgenotype")

# Turn MRN to 0, 1, 2
geno_wide_012 <- apply(geno_wide[, 2:ncol(geno_wide)], MARGIN = 2, 
                       FUN = function(x){
                         unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
                       })

geno_wide_012 <- cbind(geno_wide[,1], geno_wide_012)


# Pull Reference ID from S3

source('Credential/Credential_sdi.R')
Crop <- 'Tomato'

BaseInputBucketAddress <- "/shared/Veg_Phenotypes/Ref_Data/"
IDTableName <- paste0(BaseInputBucketAddress,Crop, "/",Crop,"_IDs.RData")
s3load(object = IDTableName, bucket = "genome-analytics-perm-space")

CropIDs_sub <- CropIDs %>% 
  dplyr::select(M.GERMPLASM.X_ID, M.GERMPLASM.PEDIGREE) %>% 
  dplyr::rename(ProgenyGermID = M.GERMPLASM.X_ID, 
         Pedigree = M.GERMPLASM.PEDIGREE) %>% 
  mutate(ProgenyGermID = as.numeric(ProgenyGermID))

geno_wide_012$ProgenyGermID <- as.numeric(geno_wide_012$ProgenyGermID)

geno_wide_ped <- CropIDs_sub %>% 
  inner_join(geno_wide_012, by = 'ProgenyGermID')

geno_wide_ped_1 <- geno_wide_ped[!duplicated(geno_wide_ped$ProgenyGermID),]


AWS.Authorization('ycao1')
s3write_using(geno_wide_ped_1, FUN = write_csv, 
              object = 'ycao1/ITG_data/ProcTomato/Geno/imputedGeno.csv', 
              bucket = bucket)
s3write_using(geno_wide_ped_1, FUN = write_csv, 
              object = '/shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/imputedGeno_wide.csv', 
              bucket = bucket)


samples_taqman <- as.numeric(unique(geno$ProgenyGermID[which(geno$ImputationStatus == 'fromLowDens')])) # 2720
samples_FP <- as.numeric(unique(geno$ProgenyGermID[which(geno$ImputationStatus == 'FPphased')])) # 582

ped_taqman <- geno_wide_ped_1[which(geno_wide_ped_1$ProgenyGermID %in% samples_taqman), c('ProgenyGermID', 'Pedigree')]
ped_FP <- geno_wide_ped_1[which(geno_wide_ped_1$ProgenyGermID %in% samples_FP), c('ProgenyGermID', 'Pedigree')]

s3write_using(ped_FP, FUN = write_csv, 
              object = 'ycao1/ITG_data/ProcTomato/Geno/ped_FP.csv', 
              bucket = bucket)
s3write_using(ped_taqman, FUN = write_csv, 
              object = 'ycao1/ITG_data/ProcTomato/Geno/ped_taqman.csv', 
              bucket = bucket)
