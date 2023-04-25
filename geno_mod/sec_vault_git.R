
load_env<-function(){
  # Dynamic way of checking required libraries in R
  list.of.packages <- c("aws.s3", "httr", "tidyverse", "asreml", "asremlPlus", "measurements", "reshape", "tidyr", 
                        "jsonlite", "dplyr", "readr", "foreach", "parallel", "doParallel","data.table","reshape2","readxl","future.apply")
  dynamic_require <- function(package){
    if(eval(parse(text=paste("require(",package,")")))){return(TRUE)}
    install.packages(package)
    return(eval(parse(text=paste("require(",package,")"))))
  }
  sapply(list.of.packages, dynamic_require)
  rename = dplyr::rename
  
  # Source functions
  source("/repos/ELBFA_GS_modules/pheno_mod/calc_trt_git.R")
  source("/repos/ELBFA_GS_modules/pheno_mod/sec_vault_git.R")
  source("/mnt/vaultCredentials.R")
  
  VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
  print(VaultToken)
  AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
  print(AWSCreds)
  Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
             "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
             "AWS_DEFAULT_REGION" = "us-east-1")
}




