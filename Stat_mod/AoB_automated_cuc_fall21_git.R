##########
# @ ENVIRONMENT : Duplicate of R 3.5 and Python 3.7 -- EVA and Asreml

install.packages('azurequest',repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))
install.packages('pedigree',repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))


# This is an example to run AoB with data pulled from H2H

# Dynamic way of checking required libraries in R
list.of.packages <- c("aws.s3", "httr", "tidyverse", "asreml", "asremlPlus", "measurements", "reshape", "tidyr", 
                      "jsonlite", "dplyr", "readr", "foreach", "parallel", "doParallel","data.table","reshape2","azurequest","pedigree")
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
crop       = 'Cucumber'

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

## example cucumber example AOB cucumber.R
## clean up with calculated traits and outlier detection of only test sets for genotyping see phenos_cucumber.R

s3load(object = '/ycao1/ITG/Cucumber/Fall/Cln_cucumber_fall12102021_new.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#load("/mnt/Cln_cucumber_fall12102021_new.Rdata")
dataset<-pheno

## UOM Conversion
UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
dataset$UOM <- UOMTable$NAME[match(dataset$UOM,UOMTable$UOM_ID)]
dataset_uom <- conv_uom(dataset) 

################################################
## load genotypes


s3load(object = '/ycao1/ITG/Cucumber/hybrids/hybrid_geno_cucumber.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
#load("/mnt/hybrid_geno_cucumber.Rdata")
g<-hybrid_wide_ped
gg = g[,2:ncol(g)]
colnames(gg)[1] = 'PEDIGREE'


## training set 2020 GROWSEASON=="2020:08"
# DH lines processing_2020_fall
s3load(object = '/ycao1/ITG/Cucumber/hybrids/Fall_2020_geno_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#load("/mnt/DH_processing/Fall_2020_geno_traincuc.Rdata")

g1<-Fall_20_train
gg1 = g1[,2:ncol(g1)]
colnames(gg1)[1] = 'PEDIGREE'

## training set 2021 GROWSEASON=="2021:08"
# DH lines processing_2021_fall
s3load(object = '/ycao1/ITG/Cucumber/hybrids/Fall_2021_geno_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

load("/mnt/DH_processing/Fall_2021_geno_traincuc.Rdata")

g2<-Fall_21_train
gg2 = g2[,2:ncol(g2)]
colnames(gg2)[1] = 'PEDIGREE'

## all genos
tst <- gg4[, 3:ncol(gg4)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
gg4[,2:ncol(gg4)] <- tst
table(is.na(gg4))

gg4=gg4%>%distinct()

#table(gg$PEDIGREE%in%dataset_uom$PEDIGREE_NAME)
#FALSE  TRUE 
#2694  1447 


################################################
# create test folder to store intermediate files for JWAS to run
#out0 = "/mnt/veg_blup/swtest"
#if(file.exists(out0)==F) {dir.create(out0, recursive = T)}

out0 = "/mnt/veg_blup/train_fall_2021_1yr/"
if(file.exists(out0)==F) {dir.create(out0, recursive = T)}

### QA/QC case 3
#remove filler

dataset_uom<-dataset_uom%>%filter(!PEDIGREE_NAME=="FILLER")

#### Run once the inbreeding code

inbred_vals(pheno_1=dataset_uom,crop='Cucumber')
#load("/mnt/5gen_inbreeding_Cucumber.Rdata")

################################################
## Calculate stage count for parents
## needed to run the master table needs to be run with the prediction code file6

#stage_count_dat <- calc_stage_count(dataset_uom)

#writing_csv <- function(dat, filename){
#  write.csv(dat, file = filename, row.names = F)
#}

#s3write_using(stage_count_dat, FUN = writing_csv,
#               object = '/Cucumber/Cucumber_STAGE_COUNTxxx.csv',
#               bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")

###################
## function SSGBLUP
###################

TRAITS<-c("EXTCO","FRLGT","AFW_C","FRNMK","SHAPE","WFSA_C","NFSA_C","FRNMKAvg_1-57","WFSA_CAvg_1-57","NFSA_CAvg_1-57","WFSA_CSum_1-57","NFSA_CSum_1-57","NETWTSum_1-57" )

### scenario 1 reduce dataset 3 years data
dataset_uom<-dataset_uom%>%filter(!GROWSEASON%in%c("2015:08","2016:08","2017:08"))

### scenario 2 training set 2021
dataset_uom<-dataset_uom%>%filter(GROWSEASON%in%c("2021:08"))

dataset_uom<-dataset_uom%>%filter(OBSRVTN_REF_CD%in%TRAITS)

#@method1 = SSGBLUP or ABLUP or GBLUP
#@seas=NULL is prediction with phenotype in model, seas="2022:08" masks the phenotype of the season
#@ gg = NULL for ABLUP


results <- data.frame()
for (i in 1:length(TRAITS)){
  trait_choosen = TRAITS[i]
  print(paste(" ################--- ", trait_choosen,"__is completed ############"))
  method1='ABLUP'
  res = getoutput(out=out0,dataset=dataset_uom,trait = trait_choosen,gg=NULL,n_gen=5,seas=NULL,method=method1)
  outbox<- paste0(out0, method1,'_deep/', trait = trait_choosen, '_deep')
  trait_GCA<-output(trait=trait_choosen,outbox=outbox,out=out0,crop=crop,method="SCA",tmp=tmp,gmethod=method1)
  results = rbind(results,trait_GCA)
  save(results,file="results_asreml.Rdata")
  #save(results,file="results_aBLUP.Rdata")
}
