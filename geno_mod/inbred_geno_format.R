library(aws.s3)
library(plyr)
library(magrittr)
library(jsonlite)
library(data.table)
library(bit64)
library(tidyr)
library(httr)

source("/repos/ELBFA_GS_modules/geno_mod/func_sql.R")
source("/repos/ELBFA_GS_modules/geno_mod/sec_vault_git.R")

load_env()
options(scipen = 999)

### training and 2020 is the same for all seasons
### 2022 is 2021 DH lines collected different for each season

### training

object1 = "Cucumber/LongDutch/Y3/SeasonSetTrainingSet_20210825_Auto/HDgenotypes.txt.gz"
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

object2 = "Cucumber/LongDutch/Y3/SeasonSetTrainingSet_20210825_Auto/CukeY3_TrainingSet_SeasonSet_MabmartFauxProgenyTable.txt"
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

object3 = 'Cucumber/LongDutch/Y3/SeasonSetTrainingSet_20210825_Auto/MarkerMetadata.txt'
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

### 2020
object1="Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/HDgenotypes.txt.gz"
bucket='veg-apd-sdi-predictiveanalytics-prod-geno-data'

object2="Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/WorkingFiles/Cucumber_LongDutch_Y3_2500Lines_ProgenyTable.txt"
bucket='veg-apd-sdi-predictiveanalytics-prod-geno-data'

object3='Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/MarkerMetadata.txt'
bucket='veg-apd-sdi-predictiveanalytics-prod-geno-data'


### 2021

### spring 2021- 2022 grown

object = "Cucumber/LongDutch/Y3/DH0_Summer2021/HDgenotypes.txt.gz"
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

object = "Cucumber/LongDutch/Y3/DH0_Summer2021/Cucumber_LongDutch_Y3_2021SummerCycleDH0_ProgenyTable.txt"
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

object = 'Cucumber/LongDutch/Y3/DH0_Summer2021/MarkerMetadata.txt'
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

### summer 2021-2022 grown

object = "Cucumber/LongDutch/Y3/DH0_Fall2021/HDgenotypes.txt.gz"
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

object = "Cucumber/LongDutch/Y3/DH0_Fall2021/Cucumber_LongDutch_Y3_2021FallCycleDH0_ProgenyTable.txt"
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

object = 'Cucumber/LongDutch/Y3/DH0_Fall2021/MarkerMetadata.txt'
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

### fall 2021

object = "Cucumber/LongDutch/Y3/DH0_Spring_complete_20210914/HDgenotypes.txt.gz"
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

object = "Cucumber/LongDutch/Y3/DH0_Spring_complete_20210914/Y3_LongDutch_2021SpringCycle_Combined_ProgenyTable.txt"
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'

object = 'Cucumber/LongDutch/Y3/DH0_Spring_complete_20210914/MarkerMetadata.txt'
bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data'


#############
### FUNCTION
plan(multisession, workers = 10) 
crop="Cucumber"
object4="s3://veg-apd-sdi-predictiveanalytcs-prod-reference-data/Cucumber/Cucumber_IDs.csv"

inbred_format=function(crop,object1,object2,object3,bucket){
  geno_DH <- s3read_using(fread, 
                          object = object1, 
                          bucket = bucket )
  
  
  meta <- s3read_using(fread, 
                       object = object2  ,
                       bucket = bucket )
  
  markerdata <- s3read_using(fread, 
                             object = object3,
                             bucket = bucket )
  
  
  # Wide formatting of the original genotypic file
  geno_wide <- geno_DH %>% 
    dplyr::select(ProgenyGermID, MRN, HDgenotype) %>% 
    dcast(ProgenyGermID ~ MRN, value.var="HDgenotype")
  
  FUN = function(x){
    unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
  }
  
  geno_wide_012<-future_apply(geno_wide[, 2:ncol(geno_wide)],MARGIN = 2, FUN=FUN)
  geno_wide<-data.frame(geno_wide)
  geno_wide_012<-data.frame(geno_wide_012)
  geno_wide_012 <- cbind(data.frame("ProgenyGermID"=geno_wide$ProgenyGermID), geno_wide_012[,1:ncol(geno_wide_012)])
  
  # Pull Reference ID from S3
  
  Crop <- crop
  obj <-get_object(object4)  
  csvcharobj <- rawToChar(obj)  
  con <- textConnection(csvcharobj)  
  CropIDs <- read.csv(file = con)
  
  CropIDs_sub <- CropIDs %>%
    dplyr::select(M.GERMPLASM.X_ID,  M.GERMPLASM.PEDIGREE) %>%
    dplyr::rename(ProgenyGermID = M.GERMPLASM.X_ID,
                  Pedigree = M.GERMPLASM.PEDIGREE) %>%
    mutate(ProgenyGermID = as.numeric(ProgenyGermID)) %>%
    unique()
  
  
  geno_wide_012$ProgenyGermID <- as.numeric(geno_wide_012$ProgenyGermID)
  
  geno_wide_ped <- CropIDs_sub %>% 
    dplyr::inner_join(geno_wide_012, by = 'ProgenyGermID')
  
  DH_data<-geno_wide_ped
}


#save(DH_2020,file="DH_2020_geno_cucumber.Rdata")