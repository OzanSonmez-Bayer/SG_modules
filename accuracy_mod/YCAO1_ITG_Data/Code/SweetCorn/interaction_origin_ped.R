# Ped(CROSS NAME):Origin interaction in the model 

library(dplyr)
library(tidyr)
packageVersion('dplyr') # 1.0.0
packageVersion('tidyr') # 1.1.0
# tidyr::pivot_longer

#install.packages('outliers', lib = "/mnt/Packages")
#install.packages('asremlPlus', lib = "/mnt/Packages")
#install.packages('measurements', lib = "/mnt/Packages")
install.packages('/mnt/packages/doBy_4.6-2.tar.gz', repos = NULL, type = 'source')

library(aws.s3)
library(plyr)
library(magrittr)
library(readxl)
library(outliers)
library(asreml)
library(asremlPlus)
library(doBy)
library(jsonlite)
library(measurements)
library(data.table)
#devtools::install_github("apache/arrow/r", lib='/mnt/Packages')
library(arrow)
library(stringr)
library(gtools)
library(httr)
library(caret)
print('Packages loaded.')


source('Credential/vaultCredentials.R')
source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')

Crop <- 'Corn'

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")


source('/repos/ITG_Data/Code/Shared/deepped.R')
# source('/repos/BLUPF90/BLUPF90/cross_validation/functStore_deep_ITG_cv.R')
# 
# loadpackage()


pmemory    = 1e9
pheno <- s3readRDS(object = 'ycao1/ITG/Corn/pheno_origin_hbc.rds',
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

ELEN_dat <- pheno %>% 
  filter(OBSRVTN_REF_CD == 'ELEN') %>% 
  mutate(TREP = paste(REP_NUMBER,TRACK_ID,sep='_')) %>% 
  mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
  mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
  mutate(pedigree = as.character(PEDIGREE_NAME)) %>% 
  mutate(BR_FIELD_ID = as.character(BR_FIELD_ID)) %>% 
  mutate(ORIGIN = as.character(ORIGIN)) %>% 
  mutate(P1_ORIGIN = as.factor(P1_ORIGIN))

bdat <- ELEN_dat %>% 
  filter(GROWSEASON %in% c('2019:04','2020:04'))

rm(ELEN_dat)
rm(pheno)
gc()

if(sum(bdat$ORIGIN == '')>0){
  bdat$ORIGIN[which(bdat$ORIGIN == '')] <- bdat$pedigree[which(bdat$ORIGIN == '')]
}

if (!'P1' %in% colnames(bdat)){
  print("parentage pedigree doesn't exit")
  bdat <- bdat %>% 
    separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
  p1 <- unlist(lapply(bdat$P1, FUN = pedSingle))
  p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  bdat$P1 <- p1
  bdat$P2[bdat$P2=='']=NA
  p2 <- unlist(lapply(bdat$P2, FUN = pedSingle))
  p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  bdat$P2 <- p2
  
}else{
  print("parentage pedigree exist")
  colnames(bdat)[which(colnames(bdat) == 'parent.female_pedigree_name')] = 'P1'
  colnames(bdat)[which(colnames(bdat) == 'parent.male_pedigree_name')] = 'P2'
  
}
if (sum(is.na(bdat$P1)) >0) {
  bdat[is.na(bdat$P1), 'P2'] <- NA
}

if (sum(is.na(bdat$P2)) >0 ){
  bdat[is.na(bdat$P2), 'P1'] <- NA
}

bdat$P1 <- as.character(bdat$P1)
bdat$P2 <- as.character(bdat$P2)

cols_chk <- c('BR_FIELD_ID', 'GROWSEASON', 'REPETITION', 'TREP')
bdat[cols_chk] <- lapply(bdat[cols_chk], factor)
levs <- sapply(bdat[cols_chk],function(x)length(unique(x)))
print(levs)
rando <- ""
if (levs['BR_FIELD_ID'] > 1){
  rp <- sapply(unlist(lapply(split(bdat,bdat$BR_FIELD_ID,drop=T),function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
  if(sum(rp, na.rm = T) >= 1){
    rando <- '~BR_FIELD_ID:TREP'
  }else{
    rando <- '~BR_FIELD_ID'
  }
}else{
  
  if (levs['TREP'] > 1){
    rando <- '~TREP'
  }
}

if (levs['GROWSEASON'] > 1){
  if (rando != ""){
    rando <- paste0(rando, '+ GROWSEASON') 
  }else{
    rando <- '~GROWSEASON'
  }
}
if (levs['REPETITION'] > 1){
  if (rando != ''){
    rando <- paste0(rando, '+ REPETITION')
  }else{
    rando <- '~REPETITION'
  }
}

if (rando != ""){
  randof <- formula(paste0(rando, '+ped(PEDIGREE_NAME)'))
}else{
  randof <- formula('~ped(PEDIGREE_NAME)')
}
print(randof)

pd <- unique(data.frame(bdat$PEDIGREE_NAME, bdat$P1, bdat$P2))
pd <- pd[!duplicated(pd[,1]),]
names(pd) <- c('Pedigree', 'P1', 'P2')
pd$Selfing <- rep(0, rep = nrow(pd))

pd.aniv=asreml.Ainverse(pd,fgen = c('Selfing',5))


met_mdl <- asreml(fixed = TRAIT_VALUE ~ 1,
                  random = ~P1_ORIGIN:ped(PEDIGREE_NAME),
                  na.method.X = 'include', trace = FALSE, 
                  ginverse = list(PEDIGREE_NAME = pd.aniv$ginv),
                  # rcov = ~at(P1_ORIGIN):units,
                  workspace = 320e+6, pworkspace = 320e+6,
                  data = bdat)

met_mdl <- update.asreml(met_mdl)

summary(met_mdl)$varcomp

pred_df <- as.data.frame(predict(met_mdl, classify = 'P1_ORIGIN:ped(PEDIGREE_NAME)', max.iter = 30,
                                 workspace = 320e+6, pworkspace = 320e+6)$predictions$pvals)

head(pred_df)
