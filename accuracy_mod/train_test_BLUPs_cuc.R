list.of.packages <- c("aws.s3", "httr", "tidyverse", "asreml", "asremlPlus", "measurements", "reshape", "tidyr", 
                      "jsonlite", "dplyr", "readr", "foreach", "parallel", "doParallel","data.table","reshape2")
dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")")))){return(TRUE)}
  install.packages(package)
  return(eval(parse(text=paste("require(",package,")"))))
}
sapply(list.of.packages, dynamic_require)
rename = dplyr::rename
#https://s3.console.aws.amazon.com/s3/object/veg-apd-sdi-predictiveanalytics-prod-geno-data?region=us-east-1&prefix=Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all_08012021.rds

######## PBLUP progeny prediction ########
#phenoTrain = phenoTrain
#phenoTest = phenotTest
#trt = trt
#crop = crop
#gg=gg
#method = method
#GCAorSCA = "SCA"


compute_accuracy <- function(method,trt, GCAorSCA = "SCA",
                             phenoTrain = NULL,
                             phenoTest = NULL,
                             gg=gg,
                             crop = 'Cucumber'
){
  if (toupper(method) == "PBLUP"){
    blup_df <- pblup_fun(phenoTrain, crop)
    respblup<-list(blup_df,phenoTest)
    save(respblup,file=paste0("/mnt/",trt,"_PBLUPresult.Rdata"))
  }
  
  if (toupper(method) == "GBLUP"){
    blup_df <- GWS_fun(trt=trt,phen=phenoTrain,method='GBLUP',crop=crop,gg=gg)
    method=method
    trt=trt
    method1=paste0("/",method,"/")
    outbox = paste0("/mnt/veg_blup/gws_cuc_test", method1, trt, '_deep')
    load("/mnt/5gen_inbreeding_comp_Cucumber.Rdata")
    trait_prediction<-output(trait=trt,outbox=outbox,crop=crop,method="none",tmp=tmp)
    res2<-list(trait_prediction,phenoTest)
    save(res2,file=paste0("/mnt/",trt,"_result.Rdata"))
  }
  
  if (toupper(method) == "SSGBLUP"){
    blup_df <- GWS_fun(trt=trt,phen=phenoTrain,method='SSGBLUP',crop=crop,gg=gg)
    method=method
    trt=trt
    method1=paste0("/",method,"/")
    outbox = paste0("/mnt/veg_blup/gws_cuc_test", method1, trt, '_deep')
    load("/mnt/5gen_inbreeding_comp_Cucumber.Rdata")
    trait_prediction<-output(trait=trt,outbox=outbox,crop=crop,method="none",tmp=tmp)
    res2<-list(trait_prediction,phenoTest)
    save(res2,file=paste0("/mnt/",trt,"_SSGBLUPresult.Rdata"))
  }
  
  
}


#function starts
makeTrainSet <- function(linesToExclude = NULL, originToExclude = NULL, pheno = pheno ){
  if (!is.null(linesToExclude)){
    pheno$TRAIT_VALUE[pheno$PEDIGREE_NAME %in% linesToExclude ] <- NA
  }
  
  if (!is.null(originToExclude)){
    pheno$TRAIT_VALUE[pheno$ORIGIN %in% originToExclude ] <- NA
  }
  return(pheno)
} 

#function ends

install.packages('https://cran.r-project.org/src/contrib/Archive/parallelly/parallelly_1.20.0.tar.gz', repos = NULL, type = 'source')
library(parallelly)

install.packages('https://cran.r-project.org/src/contrib/Archive/doBy/doBy_4.6-2.tar.gz', repos = NULL, type = 'source')
library(doBy)

source('/mnt/vaultCredentials.R')
source('/repos/YCAO1_ITG_Data/Code/Shared/BLUPFunctions.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

#s3load(object = 'GPC_Data/Corn/Corn_GPC.RData', 
#       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data')

### PHENOTYPE DATA

s3load(object = 'ycao1/ITG/Cucumber/Cln_cucumber_fall12102021_new.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#### load genos

s3load(object = 'ycao1/ITG/Cucumber/training_geno_cucumber.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
g<-hybrid_wide_ped
gg = g[,2:ncol(g)]
colnames(gg)[1] = 'PEDIGREE'

## training set 2020 GROWSEASON=="2020:08"
# DH lines processing_2020_fall
#load("/mnt/DH_processing/Fall_2020_geno_traincuc.Rdata")
#s3save(Fall_20_train, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/Fall_2020_geno_traincuc.Rdata")

s3load(object = '/ycao1/ITG/Cucumber/Fall_2020_geno_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

g1<-Fall_20_train
gg1 = g1[,2:ncol(g1)]
colnames(gg1)[1] = 'PEDIGREE'

## training set 2021 GROWSEASON=="2021:08"
# DH lines processing_2021_fall

#load("/mnt/DH_processing/Fall_2021_geno_traincuc.Rdata")
#s3save(Fall_21_train, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/Fall_2021_geno_traincuc.Rdata")
s3load(object = '/ycao1/ITG/Cucumber/Fall_2021_geno_traincuc.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

g2<-Fall_21_train
gg2 = g2[,2:ncol(g2)]
colnames(gg2)[1] = 'PEDIGREE'

## all genos
gg4=rbind(gg,gg1,gg2)
tst <- gg4[, 3:ncol(gg4)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
gg4[,2:ncol(gg4)] <- tst
table(is.na(gg4))

##### Subset phenotypes based on genotype pedigree matching

pheno=pheno%>%dplyr::filter(!GROWSEASON=="2015:08")
pheno=pheno%>%dplyr::filter(!PEDIGREE_NAME=="FILLER")
#pheno_1 <- pheno %>% 
#  mutate(year_stage = paste0(GROWSEASON, '_', EXPER_STAGE_REF_ID))

##### originToExclude in P3 P4, 2018,2019,2020
##### we exclude every observation from the year including P1

#originToExclude <- pheno_1 %>% filter(EXPER_STAGE_REF_ID %in% c('P4', 'P3','P2') & GROWSEASON %in% c("2019:08","2020:08", "2021:08")) %>% select(ORIGIN) %>% unique() %>% as.data.frame()
#hybridToExclude   <- pheno_1 %>% filter(EXPER_STAGE_REF_ID %in% c('P4', 'P3','P2') & GROWSEASON %in% c("2019:08","2020:08", "2021:08")) %>% select(PEDIGREE_NAME) %>% unique() %>% as.data.frame()
#phenoTrain_allTraits= makeTrainSet(linesToExclude = hybridToExclude$PEDIGREE_NAME, originToExclude = unlist(originToExclude$ORIGIN), pheno = pheno_1)  #remove all origins that are in the testing data
trt="AFW_C"
gg=gg4
balanced=TRUE
crop='Cucumber'
method='GBLUP'
pheno=pheno


comp_met<-function(trt,gg=gg,method,crop,pheno,balanced){
  
  ##### Subset phenotypes based on genotype pedigree matching
  pheno_1 <- pheno %>% 
    mutate(year_stage = paste0(GROWSEASON, '_', EXPER_STAGE_REF_ID))
  
  ##### originToExclude in P3 P4, 2018,2019,2020
  ##### we exclude every observation from the year including P1
  
  originToExclude <- pheno_1 %>% dplyr::filter(EXPER_STAGE_REF_ID %in% c('P4', 'P3','P2') & GROWSEASON %in% c("2019:08","2020:08", "2021:08")) %>% select(ORIGIN) %>% unique() %>% as.data.frame()
  hybridToExclude   <- pheno_1 %>% dplyr::filter(EXPER_STAGE_REF_ID %in% c('P4', 'P3','P2') & GROWSEASON %in% c("2019:08","2020:08", "2021:08")) %>% select(PEDIGREE_NAME) %>% unique() %>% as.data.frame()
  phenoTrain_allTraits= makeTrainSet(linesToExclude = hybridToExclude$PEDIGREE_NAME, originToExclude = unlist(originToExclude$ORIGIN), pheno = pheno_1)  #remove all origins that are in the testing data
  
  ## Prediction Power ###
  ### by trait
  phenoTrain = phenoTrain_allTraits %>% dplyr::filter(OBSRVTN_REF_CD == trt)
  ### subset train pheno by matching with genotype file
  if(balanced==TRUE){
    phenoTrain = phenoTrain %>%dplyr::filter(phenoTrain$PEDIGREE_NAME%in%gg$PEDIGREE)# complete dataset
  }else{
    phenoTrain = phenoTrain
  }
  ### test set phenotypes from original pheno file
  pheno_2<-pheno_1 %>% dplyr::filter(OBSRVTN_REF_CD == trt)
  phenotTest = pheno_2 %>% dplyr::filter(EXPER_STAGE_REF_ID %in% c('P4', 'P3','P2') & GROWSEASON %in% c("2019:08","2020:08", "2021:08"))

  if(balanced==TRUE){
    phenotTest = phenotTest %>%dplyr::filter(phenotTest$PEDIGREE_NAME%in%gg$PEDIGREE)
  }else{
    phenotTest = phenotTest
  }  
  
  res<-compute_accuracy( phenoTrain = phenoTrain,
                         phenoTest = phenotTest,
                         trt = trt,crop = crop,gg=gg,
                         method = method,   GCAorSCA = "SCA"
  )
  
  return(res)
}


#### output processing 
proc<-function(res,trt){
  trait_prediction<-res[[1]]
  phenoTest<-res[[2]]
  edia_pblup <- trait_prediction %>% 
    dplyr::filter(trait == trt ) 
  # %>% 
  #filter(PEDIGREE_NAME %in% unlist(lineToCompute))
  edia_raw <- phenoTest %>% 
    dplyr::filter(OBSRVTN_REF_CD == trt) %>% 
    # filter(GROWSEASON %in% growSeason) %>%
    # filter(PEDIGREE_NAME %in% unlist(lineToCompute)) %>% 
    dplyr::group_by(PEDIGREE_NAME) %>% 
    summarise(lsmean = mean(as.numeric(TRAIT_VALUE), na.rm = T))
  edia_compare <- edia_pblup %>% inner_join(edia_raw, by = 'PEDIGREE_NAME')
  edia_compare2 <<- edia_compare
  print(dim(edia_compare))
  print(cor(edia_compare$predicted.value, edia_compare$lsmean))
  res<-cor(edia_compare$predicted.value, edia_compare$lsmean)
  return(res)
}


proc_pblup<-function(res,trt){
  trait_prediction<-respblup[[1]][[2]]
  phenoTest<-respblup[[2]]
  edia_pblup <- trait_prediction %>% 
    dplyr::filter(trait == trt ) 
  # %>% 
  #filter(PEDIGREE_NAME %in% unlist(lineToCompute))
  edia_raw <- phenoTest %>% 
    dplyr::filter(OBSRVTN_REF_CD == trt) %>% 
    # filter(GROWSEASON %in% growSeason) %>%
    # filter(PEDIGREE_NAME %in% unlist(lineToCompute)) %>% 
    group_by(PEDIGREE_NAME) %>% 
    summarise(lsmean = mean(as.numeric(TRAIT_VALUE), na.rm = T))
  edia_compare <- edia_pblup %>% 
    inner_join(edia_raw, by = 'PEDIGREE_NAME')
  edia_compare2 <<- edia_compare
  print(dim(edia_compare))
  print(cor(edia_compare$predicted.value, edia_compare$lsmean))
  res<-cor(edia_compare$predicted.value, edia_compare$lsmean)
  return(res)
}



###################
#### GBLUP
###################


trait<-c("EXTCO","FRLGT","AFW_C","FRNMK","SHAPE","WFSA_C","NFSA_C")
trait<-c("FRNMKAvg_1-57","WFSA_CAvg_1-57","NFSA_CAvg_1-57")
trait<-c("EXTCO")


library(parallel)
mclapply(trait,function(x) comp_met(x,gg=gg4,method='SSGBLUP',crop='Cucumber',pheno=pheno,balanced=TRUE), mc.cores = 5,mc.preschedule = FALSE)
mclapply(trait,function(x) comp_met(x,gg=gg4,method="PBLUP",crop='Cucumber',pheno=pheno,balanced=TRUE), mc.cores = 5,mc.preschedule = FALSE)
mclapply(trait,function(x) comp_met(x,gg=gg4,method="GBLUP",crop='Cucumber',pheno=pheno,balanced=TRUE), mc.cores = 3,mc.preschedule = FALSE)


mclapply(trait,function(x) comp_met(x,gg=gg4,method='SSGBLUP',crop='Cucumber',pheno=pheno,balanced=FALSE), mc.cores = 3,mc.preschedule = FALSE)

mclapply(trait,function(x) comp_met(x,gg=gg4,method="PBLUP",crop='Cucumber',pheno=pheno,balanced=TRUE), mc.cores = 5,mc.preschedule = FALSE)

mclapply(trait,function(x) comp_met(x,gg=gg4,method="GBLUP",crop='Cucumber',pheno=pheno,balanced=TRUE), mc.cores = 3,mc.preschedule = FALSE)

