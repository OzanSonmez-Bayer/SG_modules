# This is an example to run AoB with data pulled from H2H
#trt=trt
#phen=phenoTrain
#method='GBLUP'
#crop=crop
#gg=gg

GWS_fun<-function(phen,method,crop,trt,gg){
  # Dynamic way of checking required libraries in R
  list.of.packages <- c("aws.s3", "httr", "tidyverse", "asreml", "asremlPlus", "measurements", "reshape", "tidyr", 
                        "jsonlite", "dplyr", "readr", "foreach", "parallel", "doParallel","data.table","reshape2")
  dynamic_require <- function(package){
    if(eval(parse(text=paste("require(",package,")")))){return(TRUE)}
    install.packages(package)
    return(eval(parse(text=paste("require(",package,")"))))
  }
  sapply(list.of.packages, dynamic_require)
  rename = dplyr::rename
  
  # Source functions
  source("/mnt/deep_ped_functions_ODB.R")
  source('/mnt/vaultCredentials.R')
  source('/mnt/ValutCreds_ancestry.R')
  source('/mnt/getbdat_AOB_V3.R')
  source('/mnt/JWAS_Model.R')
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
  
  ################################################
  #### read cropIDs
  crop       = crop
  IDTableName <- paste0(crop, '/', crop, '_IDs.RData')
  print(IDTableName)
  s3load(object = IDTableName, bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")
  cropids    = CropIDs[, c('M.GERMPLASM.X_ID','M.GERMPLASM.PEDIGREE')]
  cropids    = cropids[duplicated(cropids)==F, ]
  
  ################################################
  #### read phenos
  ## example corn
  #dataset <- s3read_using(readRDS, bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace', 
  #                        object = 'ycao1/ITG/Corn/pheno.rds')
  dataset <- phen ## validation phenotype
  
  ## UOM Conversion
  UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
  dataset$UOM <- UOMTable$NAME[match(dataset$UOM,UOMTable$UOM_ID)]
  dataset_uom <- conv_uom(dataset) 
  
  ################################################
  ## load genotypes
  #s3load(object = 'ycao1/ITG/Cucumber/training_geno_cucumber.Rdata', 
  #       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
  #g<-hybrid_wide_ped
  #gg = g[,2:ncol(g)]
  #colnames(gg)[1] = 'PEDIGREE'
  
  
  ################################################
  # create test folder to store intermediate files for JWAS to run
  out0 = "/mnt/veg_blup/gws_cuc_test"
  if(file.exists(out0)==F) {dir.create(out0, recursive = T)}
  
  
  #### Run once the 5gen_inbreeding_SweetCorn.Rdata 
  inbred_vals(pheno_1=dataset,crop=crop)
  load(paste0("/mnt/5gen_inbreeding_comp_",crop,".Rdata"))
  
  
  ###################
  ## function SSGBLUP
  ###################
  
  ## results in long format
  
  res = getoutput( method=method,out=out0,dataset=dataset_uom,trait = trt,gg=gg,n_gen=5,seas=NULL,crop=crop)
  #outbox<-paste0(out0, '/SSGBLUP_deep/', trait=trait_choosen, '_deep')
}

