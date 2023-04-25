# Bayesian based GBLUPs

library(aws.s3)
install.packages('BGLR')
library(BGLR)
library(httr)

source('/mnt/vaultCredentials.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

########### Import Data ##########

geno <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds', 
                  bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

pheno <- s3readRDS(object = 'ycao1/ITG/Corn/pheno_1.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

