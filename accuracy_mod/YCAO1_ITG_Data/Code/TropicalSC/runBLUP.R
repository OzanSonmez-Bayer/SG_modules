# This is an example to run AoB with data pulled from H2H

library(aws.s3)
library(httr)
library(tidyverse)
install.packages("azurequest", repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))
library(azurequest)
library(asreml)
library(asremlPlus)
install.packages('packages/doBy_4.6-2.tar.gz', type = 'source', repos = NULL)
library(doBy)
library(measurements)
library(parallel)
library(aws.s3)
library(httr)


# source('VaultCreds.R')
# source('ValutCreds_ancestry.R')

source('VaultCredentials.R')
source('credentials.R')
source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')
source('/repos/ITG_Data/Code/Shared/misFunctions.R')


VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")


Sys.setenv("CLIENT_ID" = CLIENT_ID,
           "CLIENT_SECRET" = CLIENT_SECRET)


# Define paramters 

args <- commandArgs(trailingOnly = TRUE)
Crop <- args[1]
# PMC <- args[2] # Comma separated 
traits <- args[2] # Comma separated
years <- args[3] # Comma separated
mdl <- args[4]
stages <- args[5]
output <- args[6]
savefile <- args[7]

if(file.exists(output)==F)  {dir.create(output, recursive = T)}


############### Step 1: Pull Phenotypic Data #################

pheno <- s3readRDS(object = 'ycao1/AoB/SweetCorn/Tropical/rawdata/pheno_2018_2021_clean_08112021.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

############# Step 2: Calculate Stage Count ############
stage_count_dat <- calc_stage_count(pheno)
s3write_using(stage_count_dat, FUN = writing_csv, object = paste0(output,'stage_count_dat.csv'),
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


############## Step 3: Data QA/QC #####################

# source('/repos/BLUP_VEG/veg_blup/miscFunctions.R')
# source('/repos/BLUP_VEG/veg_blup/deep_pedigree_functions.R') 
source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')
source('/repos/ITG_Data/Code/Shared/misFunctions.R')

# Note function from master branch is the up-to-date. 
# https://github.platforms.engineering/VegetablePipelineSolutions/BLUP_VEG_Master/tree/master/veg_blup

# Check and convert UOM 

# Note: column names need to be changed to H2H column names

UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
# pheno$UOM <- UOMTable$NAME[match(pheno$UOM,UOMTable$UOM_ID)]
pheno <- conv_uom(pheno)


############## Step 4: Run BLUPs ######################

# source('/repos/BLUP_VEG/veg_blup/iblup.r')
source('/repos/ITG_Data/Code/Shared/deepped.R')
source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')

# Note: column names need to be changed to H2H column names

# Deep pedigree based BLUP

# blupResults <- Ablup_fun(dataset = AllData_1, year = NA, n_gen = 5, CropIDs = CropIDs, max_memory = 320e+6)

stage <- trimws(unlist(strsplit(stages, split = ',')))
season <- trimws(unlist(strsplit(years, split = ',')))

pheno <- pheno %>% 
  # filter(EXPER_STAGE_REF_ID %in% stage) %>% 
  filter(GROWSEASON %in% season)

# traits <- c('EDIA', 'EHT', 'ELEN', 'ERSC', 'EVG', 'HSC', 'KRN', 
#             'LDSR_SPIRKU', 'LDSR_PHSPMA', 'LDSR_SETOTU', 'P50D', 'PCTRC', 
#             'PHT', 'QUAL', 'RTLP', 'S50D', 'SC_FR', 'SC_HE', 'SC_TF', 'SCCSA', 'SCTPA')

# traits <- c('LDSR_SETOTU', 'P50D', 'PCTRC', 
#             'PHT', 'QUAL', 'RTLP', 'S50D', 'SC_FR', 'SC_HE', 'SC_TF', 'SCCSA', 'SCTPA')

for (tt in traits){
  print(tt)
  temp <- pheno %>% 
    filter(OBSRVTN_REF_CD == tt)
  if (nrow(temp) < 50){
    next
  }
  if (mdl == 'ABLUP'){
    blupResults <- Ablup_fun(dataset = temp,  n_gen = 5, CropIDs = CropIDs, crop = Crop)
    # write.csv(blupResults[[1]], file = paste0(output, '/', tt, '_GCA.csv'), row.names = F)
    # write.csv(blupResults[[2]], file = paste0(output, '/', tt, '_SCA.csv'), row.names = F)
    
    s3write_using(blupResults[[1]], write.csv, 
                  object = paste0(output, savefile, '/', tt, '_GCA.csv'), 
                  bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
    
    s3write_using(blupResults[[2]], write.csv, 
                  object = paste0(output, savefile, '/', tt, '_SCA.csv'), 
                  bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
    
    rm(temp)
    rm(blupResults)
    gc()
    
  }

  if(mdl == 'BLUE'){
    print(length(unique(temp$CROSS_NAME)))
    
    NumberOfCrossNames <- length(unique(temp$CROSS_NAME))
    
    CrossesPerRun <- ceiling(NumberOfCrossNames/10)
    
    if(CrossesPerRun > 1500) { 
      CrossesPerRun = 1500 
      }
    if (CrossesPerRun < 500) { 
      CrossesPerRun = NumberOfCrossNames 
      }
    
    
    ListCrossNames <- unique(temp$CROSS_NAME)
    #ListOfList <- split(x=ListCrossNames, sort(ListCrossNames%%CrossesPerRun))
    ListOfList <- split(ListCrossNames,as.numeric(gl(length(ListCrossNames),CrossesPerRun,length(ListCrossNames))))
    
    print(length(ListOfList))
    print(CrossesPerRun)
    
    if(length(grep("Windows",Sys.getenv("OS")))>0) {
      NumberOfCores = 1
    } else {
      NumberOfCores = detectCores()/2
      
      print(paste0("Number of cores detected ",detectCores()))
      
    }
    
    if(NumberOfCores < 1) { NumberOfCores = 1}
    
    ListOf_PredBlue <- mclapply(ListOfList, function(x) RunBLUE(temp,x), mc.cores = NumberOfCores, mc.silent = FALSE)
    
    pred_blu <- do.call(rbind,ListOf_PredBlue)
    pred_blu$trait <- tt
    
    s3write_using(pred_blu, writing_csv, 
                  object = paste0(output, savefile, '/', tt, '_BLUE.csv'), 
                  bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
    
    rm(pred_blu)
    gc()
    
  }

}
# blupResults <- Ablup_fun(dataset = pheno,  n_gen = 5, CropIDs = CropIDs, crop = Crop)
# write.csv(blupResults[[1]], file = paste0(output, '/_GCA.csv'), row.names = F)
# write.csv(blupResults[[2]], file = paste0(output, '/_SCA.csv'), row.names = F)



############# Step 5: Combine All the BLUPs ###############
# res_bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-workspace'
# path <- 'ycao1/OriginSelection/Tomato/CCT_HE/ABLUP/'
# group <- 'PCMGroup4'
# 
# my_bucket <- get_bucket(bucket = res_bucket, prefix = paste0(path, group, '/'))
# 
# GCA_df <- data.frame()
# SCA_df <- data.frame()
# 
# for (list_file in my_bucket){
#   
#   print(list_file$Key)
#   
#   # generate a GCA master file 
#   
#   if (grepl('_GCA.csv', list_file$Key)){
#     temp <- s3read_using(FUN = read.csv, 
#                          object = list_file$Key, 
#                          bucket = res_bucket)
#     
#     head(temp)
#     
#     temp <- temp %>% 
#       select(PEDIGREE_NAME, predicted.value, standard.error, BLUP, BLUP.se, reliability, n, N, h2, trait)
#     temp <- temp[!duplicated(temp$PEDIGREE_NAME), ]
#     temp <- temp[!is.na(temp$PEDIGREE_NAME),]
#     
#     colnames(temp)[which(!colnames(temp) %in% c('PEDIGREE_NAME', 'trait'))] <- paste(unique(temp$trait), 
#                                                                                      colnames(temp)[which(!colnames(temp) %in%c('PEDIGREE_NAME', 'trait'))],
#                                                                                      sep = '_')
#     
#     temp <- temp %>% 
#       select(-trait)
#     
#     if(nrow(GCA_df) == 0){
#       GCA_df <- temp
#     }else{
#       GCA_df <- GCA_df %>% 
#         full_join(temp, by = 'PEDIGREE_NAME')
#     }
#   }
#   
#   # generate SCA master file
#   
#   if (grepl('_SCA.csv', list_file$Key)){
#     temp <- s3read_using(FUN = read.csv, 
#                          object = list_file$Key, 
#                          bucket = res_bucket)
#     
#     head(temp)
#     
#     temp <- temp %>% 
#       select(PEDIGREE_NAME, predicted.value, standard.error, BLUP, BLUP.se, reliability, n, N, h2, trait)
#     temp <- temp[!duplicated(temp$PEDIGREE_NAME), ]
#     temp <- temp[!is.na(temp$PEDIGREE_NAME),]
#     
#     colnames(temp)[which(!colnames(temp) %in% c('PEDIGREE_NAME', 'trait'))] <- paste(unique(temp$trait), 
#                                                                                      colnames(temp)[which(!colnames(temp) %in%c('PEDIGREE_NAME', 'trait'))],
#                                                                                      sep = '_')
#     
#     temp <- temp %>% 
#       select(-trait)
#     
#     if(nrow(SCA_df) == 0){
#       SCA_df <- temp
#     }else{
#       SCA_df <- SCA_df %>% 
#         full_join(temp, by = 'PEDIGREE_NAME')
#     }
#   }
# }
# 
# 
# # save to S3 
# 
# 
#   
# s3write_using(GCA_df, FUN = write.csv, 
#                 object = paste0(path, 'masterfile/', group, '_GCA_master.csv'),
#                 bucket = res_bucket)
# s3write_using(SCA_df, FUN = write.csv, 
#               object = paste0(path, 'masterfile/', group, '_SCA_master.csv'),
#               bucket = res_bucket)
# 

############# Step 6: Master Table ##################

# Add GCV/marker information

# Need a code review from Sebastian and Xing for the generation of master table