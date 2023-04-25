master_table<-function(){
library("aws.s3")
library(data.table)
source('/mnt/vaultCredentials.R')
source('/mnt/ValutCreds_ancestry.R')

options(scipen = 999)

#### Vault credentials
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

Sys.setenv("CLIENT_ID" = CLIENT_ID,
           "CLIENT_SECRET" = CLIENT_SECRET)


master_python<-function(file_name,file1,file2,file3,file4,file5,file6,S3){
  require(reticulate)
  use_python("/usr/local/lib/python3.7.7/")
  source_python("master_table_forR.py")
  
  if(S3==FALSE & file.exists("file_name")==FALSE){
    blup<-read.csv(file1)
    blup = create_blup_wide(blup, 'trait', 'PEDIGREE_NAME')
  }else{
    blup<-file_blup(file_name)
  }
  
  id_table<-combine_table(file2)
  master_file = combine_blup_id(blup, id_table)
  gpc_df = combine_fileobs(file3)
  master_file = combine_blup_obs(master_file, gpc_df)    
  marker_data = combine_filemarkerData(file4)
  master_file = combine_blup_markerData(master_file, marker_data)
  parentalTable = combine_fileparentalTable(file5)
  master_file = combine_blup_parental(master_file, parentalTable)
  stageData = combine_filestageData(file6)
  master_file = combine_blup_stage(master_file, stageData)
  return(master_file)
}


#setwd("/repos/ELBFA_GS_modules/diallel_mod/")
#file_name = 'veg-apd-sdi-predictiveanalytcs-prod-workspace/ycao1/OriginSelection/Corn/ThreeYear_BLUP/GCA/Temperate_ThreeYear_2020.csv' # from AOB function
file1 ='/repos/ELBFA_GS_modules/diallel_mod/GCA_master_test_sub.csv'# from AOB function
file2='veg-apd-sdi-predictiveanalytcs-prod-reference-data/Cucumber/Cucumber_IDs.csv' #available
file3='veg-apd-sdi-predictiveanalytcs-prod-pheno-data/GPC_Data/Cucumber/Cucumber_GPC.csv'
file4 ='veg-apd-sdi-predictiveanalytcs-prod-pheno-data/MarkerPredicted_Data/Cucumber/Cucumber_MarkerData.csv'
file5 = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data/Cucumber/Cucumber_Parentals.csv'
file6 = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data/Cucumber/Cucumber_STAGE_COUNT.csv' # from AOB function

# function to transform from long to wide format and to return a master table for origin prediction

master_test<-master_python(file_name,file1,file2,file3,file4,file5,file6,S3="FALSE")
fwrite(master_test, file ="/repos/ELBFA_GS_modules/diallel_mod/master_wide_test_cuc.csv")

# write to S3 and run origin prediction
test_cuc<-read.csv(file="/repos/ELBFA_GS_modules/diallel_mod/master_wide_test_cuc.csv")

writing_csv <- function(dat, filename){
  write.csv(dat, file = filename, row.names = F)
}

s3write_using(test_cuc, FUN = writing_csv,
              object = '/ycao1/OriginSelection/Test_master_cuc_git.csv',
              bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace")


print("EOF master table")
}

master_table()


#library(reticulate)
#use_python("/usr/local/anaconda/lib/python3.7/")
#source_python("master_table_forR.py")
#file_name = 'veg-apd-sdi-predictiveanalytcs-prod-workspace/ycao1/OriginSelection/Corn/ThreeYear_BLUP/GCA/Temperate_ThreeYear_2020.csv'
#blup<-file_blup(file_name)
#file2='veg-apd-sdi-predictiveanalytcs-prod-reference-data/Corn/Corn_IDs.csv'
#id_table<-combine_table(file2)
#master_file = combine_blup_id(blup, id_table)
#file3='veg-apd-sdi-predictiveanalytcs-prod-pheno-data/GPC_Data/Corn/Corn_GPC.csv'
#gpc_df = combine_fileobs(file3)
#master_file = combine_blup_obs(master_file, gpc_df)    
#file4 = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data/Corn/FTS_MarkerData_SweetCorn_2019.csv'
#marker_data = combine_filemarkerData(file4)
#master_file = combine_blup_markerData(master_file, marker_data)
#file5 = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data/Corn/Corn_Parentals.csv'
#parentalTable = combine_fileparentalTable(file5)
#master_file = combine_blup_parental(master_file, parentalTable)
#file6 = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data/Corn/Corn_FRESH_STAGE_COUNT.csv'
#stageData = combine_filestageData(file6)
#master_file = combine_blup_stage(master_file, stageData)
