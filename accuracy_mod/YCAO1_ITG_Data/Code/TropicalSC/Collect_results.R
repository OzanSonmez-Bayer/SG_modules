# collect resutls 

library(tidyverse)
library(aws.s3)
library(httr)


source('VaultCredentials.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

pheno <- s3readRDS(object = 'ycao1/AoB/SweetCorn/Tropical/rawdata/pheno_2018_2021_clean_08112021.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


res_bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-workspace'
path <- 'ycao1/AoB/SweetCorn/Tropical/BLUE/'
group <- 'ThreeYear'

my_bucket <- get_bucket(bucket = res_bucket, prefix = paste0(path, group, '/'))

GCA_df <- data.frame()
SCA_df <- data.frame()
BLUE_df <- data.frame()

for (list_file in my_bucket){

  print(list_file$Key)

  # generate a GCA master file

  if (grepl('_GCA.csv', list_file$Key)){
    temp <- s3read_using(FUN = read.csv,
                         object = list_file$Key,
                         bucket = res_bucket)

    head(temp)

    temp <- temp %>%
      select(PEDIGREE_NAME, predicted.value, standard.error, BLUP, BLUP.se, reliability, n, N, h2, trait)
    temp <- temp[!duplicated(temp$PEDIGREE_NAME), ]
    temp <- temp[!is.na(temp$PEDIGREE_NAME),]

    colnames(temp)[which(!colnames(temp) %in% c('PEDIGREE_NAME', 'trait'))] <- paste(unique(temp$trait),
                                                                                     colnames(temp)[which(!colnames(temp) %in%c('PEDIGREE_NAME', 'trait'))],
                                                                                     sep = '_')

    # colnames(temp)[which(!colnames(temp) %in% c('PEDIGREE_NAME', 'trait'))] <- paste(
    #                                                                                  colnames(temp)[which(!colnames(temp) %in%c('PEDIGREE_NAME', 'trait'))],
    #                                                                                  unique(temp$trait),
    #                                                                                  sep = '_')

    temp <- temp %>%
      select(-trait)

    if(nrow(GCA_df) == 0){
      GCA_df <- temp
    }else{
      GCA_df <- GCA_df %>%
        full_join(temp, by = 'PEDIGREE_NAME')
    }
  }

  # generate SCA master file

  if (grepl('_SCA.csv', list_file$Key)){
    temp <- s3read_using(FUN = read.csv,
                         object = list_file$Key,
                         bucket = res_bucket)

    head(temp)

    temp <- temp %>%
      select(PEDIGREE_NAME, predicted.value, standard.error, BLUP, BLUP.se, reliability, n, N, h2, trait)
    temp <- temp[!duplicated(temp$PEDIGREE_NAME), ]
    temp <- temp[!is.na(temp$PEDIGREE_NAME),]
# 
#     colnames(temp)[which(!colnames(temp) %in% c('PEDIGREE_NAME', 'trait'))] <- paste(
#                                                                                      colnames(temp)[which(!colnames(temp) %in%c('PEDIGREE_NAME', 'trait'))],
#                                                                                      unique(temp$trait),
#                                                                                      sep = '_')

    
    colnames(temp)[which(!colnames(temp) %in% c('PEDIGREE_NAME', 'trait'))] <- paste(
      
      unique(temp$trait),
      
      colnames(temp)[which(!colnames(temp) %in%c('PEDIGREE_NAME', 'trait'))],
      
      sep = '_')
    
    temp <- temp %>%
      select(-trait)

    if(nrow(SCA_df) == 0){
      SCA_df <- temp
    }else{
      SCA_df <- SCA_df %>%
        full_join(temp, by = 'PEDIGREE_NAME')
    }
  }
  
  if (grepl('_BLUE.csv', list_file$Key)){
    temp <- s3read_using(FUN = read.csv,
                         object = list_file$Key,
                         bucket = res_bucket)
    
    head(temp)
    
    temp <- temp %>%
      select(CROSS_NAME, predicted.value, standard.error, trait) %>% 
      mutate(CROSS_NAME = as.character(CROSS_NAME))
    temp <- temp[!duplicated(temp$CROSS_NAME), ]
    temp <- temp[!is.na(temp$CROSS_NAME),]
    
    trt <- strsplit(gsub('_BLUE.csv', '', list_file$Key), split = '/')[[1]]
    trt <- trt[length(trt)]
    
    pheno_sub <- pheno %>% 
      filter(OBSRVTN_REF_CD == trt)
    
    n_obs <- summaryBy(TRAIT_VALUE ~ CROSS_NAME, data = pheno_sub, FUN = function(x) {length(x)})
    colnames(n_obs)[2] <- 'n'
    
    temp <- temp %>% 
      left_join(n_obs, by = 'CROSS_NAME') %>% 
      mutate(N = nrow(pheno_sub))
  
    
    # colnames(temp)[which(!colnames(temp) %in% c('PEDIGREE_NAME', 'trait'))] <- paste(unique(temp$trait),
    #                                                                                  colnames(temp)[which(!colnames(temp) %in%c('PEDIGREE_NAME', 'trait'))],
    #                                                                                  sep = '_')
    # 
    colnames(temp)[which(!colnames(temp) %in% c('CROSS_NAME', 'trait'))] <- paste(
      
      unique(temp$trait),
      
      colnames(temp)[which(!colnames(temp) %in%c('CROSS_NAME', 'trait'))],
      
      sep = '_')
    
    temp <- temp %>%
      select(-trait)
    
    if(nrow(BLUE_df) == 0){
      BLUE_df <- temp
    }else{
      BLUE_df <- BLUE_df %>%
        full_join(temp, by = 'CROSS_NAME')
    }
  }
}



# Add HBC

s3load(object = 'GPC_Data/Corn/Corn_GPC.RData', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data')

SCA_df <- SCA_df %>% 
  left_join(GCP_TableView_noProg[, c('Pedigree', 'HBC')], by = c('PEDIGREE_NAME' = 'Pedigree'))

GCA_df <- GCA_df %>% 
  left_join(GCP_TableView_noProg[, c('Pedigree', 'HBC')], by = c('PEDIGREE_NAME' = 'Pedigree'))

BLUE_df <- BLUE_df %>% 
  left_join(GCP_TableView_noProg[, c('CrossName', 'HBC')], by = c('CROSS_NAME' = 'CrossName'))

# Add IS_CHK to SCA_df

pheno <- s3readRDS(object = 'ycao1/AoB/SweetCorn/Tropical/rawdata/pheno_2018_2021_clean.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


chk_ls <- unique(pheno$PEDIGREE_NAME[pheno$IS_CHK == 'T'])
chk_cn_ls <- unique(pheno$CROSS_NAME[pheno$IS_CHK == 'T'])

SCA_df <- SCA_df %>% 
  mutate(IS_CHK = NA, 
         IS_CHK = if_else(PEDIGREE_NAME %in% chk_ls, T, F))

BLUE_df <- BLUE_df %>% 
  mutate(IS_CHK = NA, 
         IS_CHK = if_else(CROSS_NAME %in% chk_cn_ls, T, F))

# Save results into S3

s3write_using(SCA_df, writing_csv, 
              object = paste0(path, 'ABLUP_SCA_df_threeyear.csv') ,
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

s3write_using(GCA_df, writing_csv, 
              object = paste0(path, 'ABLUP_GCA_df_threeyear.csv') ,
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

s3write_using(BLUE_df, writing_csv, 
              object = paste0(path, 'BLUE_df_threeyear.csv') ,
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

