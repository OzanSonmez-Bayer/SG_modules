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
source("/repos/ELBFA_GS_modules/Stat_mod/deep_ped_functions_ODB.R")
source('/mnt/vaultCredentials.R')
source('/mnt/ValutCreds_ancestry.R')
source('/repos/ELBFA_GS_modules/formatf90_mod/getbdat_AOB_V3_dataformat_f90_git.R')
source('/repos/ELBFA_GS_modules/formatf90_mod/JWAS_Model_dataformat_gca_f90_git.R')
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
## example cucumber example AOB cucumber.R
## clean up with calculated traits and outlier detection of only test sets for genotyping see phenos_cucumber.R

s3load(object = '/ycao1/ITG/Cucumber/Summer/Cln_cucumber_summer18022022_new.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#load("/mnt/Cln_cucumber_summer18022022_new.Rdata")
dataset<-pheno

## UOM Conversion
UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
dataset$UOM <- UOMTable$NAME[match(dataset$UOM,UOMTable$UOM_ID)]
dataset_uom <- conv_uom(dataset) 

################################################
## load genotypes

## training set 2015-2019 might need reprocessing from the genomics team if names had been changed
#load("/mnt/hybrid_geno_cucumber.Rdata")

s3load(object = '/ycao1/ITG/Cucumber/hybrids/hybrid_geno_cucumber.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

g<-hybrid_wide_ped
gg = g[,2:ncol(g)]
colnames(gg)[1] = 'PEDIGREE'

## training set 2020
s3load(object = '/ycao1/ITG/Cucumber/hybrids/summer_20_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#load("/mnt/summer_2022_DHxTsr/summer_20_traincuc.Rdata")
g1<-summer_20_train
gg1 = g1[,2:ncol(g1)]
colnames(gg1)[1] = 'PEDIGREE'

## training set 2021
s3load(object = '/ycao1/ITG/Cucumber/hybrids/summer_21_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#load("/mnt/summer_2022_DHxTsr/summer_21_traincuc.Rdata")
g2<-summer_21_train
gg2 = g2[,2:ncol(g2)]
colnames(gg2)[1] = 'PEDIGREE'

## training set 2022
s3load(object = '/ycao1/ITG/Cucumber/hybrids/summer_22_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
#load("/mnt/summer_2022_DHxTsr/summer_22_traincuc.Rdata")
g3<-summer_22_train
gg3 = g3[,2:ncol(g3)]
colnames(gg3)[1] = 'PEDIGREE'


#table(gg$PEDIGREE%in%dataset_uom$PEDIGREE_NAME)
#FALSE  TRUE 
#2694  1447 

###################
## inbreds
#summer_DH,DH_2020,inbred_geno
#load("/mnt/DH_processing/DH_2020_geno_cucumber.Rdata")#2020 DH2020
#load("/mnt/inbred_geno_cucumber.Rdata")#training inbred_geno
#load("/mnt/DH_processing/summer_DH_geno_cucumber.Rdata")#2021 summer_DH


s3load(object = '/ycao1/ITG/Cucumber/inbreds/DH_2020_geno_cucumber.Rdata', #2020 DH2020
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
s3load(object = '/ycao1/ITG/Cucumber/inbreds/inbred_geno_cucumber.Rdata', #training inbred_geno
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
s3load(object = '/ycao1/ITG/Cucumber/inbreds/summer_DH_geno_cucumber.Rdata', #2021 summer_DH
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

inbred_geno=dplyr::rename(inbred_geno,c("ProgenyGermID"="GermID") )
inb=rbind(summer_DH,DH_2020,inbred_geno)
inb=rename(inb,c('PEDIGREE'='Pedigree'))
inb=inb[,2:ncol(inb)]

## all genos
gg4=rbind(gg,gg1,gg2,gg3,inb)

tst <- gg4[, 3:ncol(gg4)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
gg4[,2:ncol(gg4)] <- tst
table(is.na(gg4))

################################################
# create test folder to store intermediate files for JWAS to run

out0 = "/mnt/veg_blup/traintest_output/"
if(file.exists(out0)==F) {dir.create(out0, recursive = T)}

### QA/QC case 3
#remove filler
dataset_uom<-dataset_uom%>%filter(!PEDIGREE_NAME=="FILLER")

###################
## function SSGBLUP
###################
# @ This function will create the input data needed to run the Single Step run using the python wrapper
# @ it has not been tested to run other blupf90 methods like GBLUP or ABLUP  
# loop function
#check traits
unique(dataset_uom$OBSRVTN_REF_CD)
TRAITS<-c("EXTCO","FRLGT","AFW_C","FRNMK","SHAPE","WFSA_C","NFSA_C","FRNMKAvg_1-76","WFSA_CAvg_1-76","NFSA_CAvg_1-76","NETWTAvg_1-76", "NFSA_CSum_1-76","WFSA_CSum_1-76","NETWTSum_1-76")
dataset_uom<-dataset_uom%>%filter(OBSRVTN_REF_CD%in%TRAITS)

### scenario 1 training set 2018-2021
dataset_uom<-dataset_uom%>%filter(!GROWSEASON%in%c("2015:05","2016:05","2017:05"))

### scenario 2 training set 1 yr BLUPs
#dataset_uom<-dataset_uom%>%filter(GROWSEASON%in%c("2021:05"))


#### 
# This code will run all the traits in TRAIT vector and output files to input in blupf90 single step
# season= NULL will not mask the season's phenotype otherwise TRAIT_VALUE=-999 which is in the options value assigned in the blupf90 parameter file

for (i in 1:length(TRAITS)){
  trait_choosen = TRAITS[i]
  getoutput_f90(out=out0,dataset=dataset_uom,trait = trait_choosen, gg=gg4,n_gen=5,seas=NULL)#
  print(paste(" ################--- ", trait_choosen,"__is completed ############"))
}



