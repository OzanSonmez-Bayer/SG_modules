####
source("/repos/ELBFA_GS_modules/geno_mod/sec_vault_git.R")
source("/repos/ELBFA_GS_modules/Stat_mod/getbdat_AOB_V3_git.R")
load_env()

#how results folder fall 2021 was uploaded to 
#lapply(dir("/mnt/veg_blup/train_fall_2021/SSGBLUP_deep",full.names = TRUE, recursive = TRUE), function(filename) {
#  put_object(file = filename, object = filename, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace")
#})

#### collect results in wide format

path="/ycao1/ITG/Cucumber/GWS_results/veg_blup/train_fall_2021/SSGBLUP_deep/"
trait_choosen<-c("AFW_C","EXTCO","FRLGT","FRNMK","SHAPE","WFSA_C","NFSA_C","FRNMKAvg_1-57","WFSA_CAvg_1-57","NFSA_CAvg_1-57","NFSA_CSum_1-57", "WFSA_CSum_1-57", "NETWTSum_1-57" )
res2=data.frame()

for (i in 1:length(trait_choosen)){
  s3load(object = paste0(path,trait_choosen[i],'_deep','/results_asreml.Rdata'), 
         bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
  res1<-results%>%filter(trait==trait_choosen[i])
  colnames(res1)<-paste(colnames(res1[,1:ncol(res1)]),trait_choosen[i],sep="_")
  colnames(res1)[1]<-c("PEDIGREE_NAME")
  if (i==1){
    res2=rbind(res2,res1)
  }else{
    res2<-left_join(res2,res1)
  }
  
}

################################
#### function reference Ozan ###

diallel<-function(trait,results,tester,inbr){
  dial_df<-expand.grid(tester, inbr)
  tom_trait<-results%>%dplyr::select(PEDIGREE_NAME,paste(c("predicted.value","standard.error","BLUP.se", "BLUP", "reliability", "trait","N", "n"  , "h2", "mdl" ),trait,sep="_")
  )
  colnames(dial_df) = c("P1", "P2")
  dial_df <- dial_df %>% left_join(tom_trait %>% select(PEDIGREE_NAME,matches( "predicted.value")), by = c("P1" = "PEDIGREE_NAME")) %>% 
    dplyr::rename(`P1 GEBV` = paste("predicted.value",trait,sep="_"))
  dial_df <- dial_df %>% left_join(tom_trait %>% select(PEDIGREE_NAME, matches( "predicted.value")), by = c("P2" = "PEDIGREE_NAME")) %>% dplyr::rename(`P2 GEBV` = paste("predicted.value",trait,sep="_"))
  nam<-paste0("MidParent_GEBV_",trait)
  dial_df <- dial_df %>% dplyr::mutate(!!nam := (`P1 GEBV` + `P2 GEBV`)/2)
  colnames(dial_df)[3:4] = paste0(trait, "_", colnames(dial_df)[3:4])
  return(dial_df) 
}

######################
#### reference file Fall dataset
s3load(object = '/ycao1/ITG/Cucumber/Fall/Cln_cucumber_fall12102021_new.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
dataset<-pheno

## UOM Conversion
UOMTable <- read.csv("/mnt/MIDAS_UOM_Table.csv", colClasses = rep(times=5,x="character"))
dataset$UOM <- UOMTable$NAME[match(dataset$UOM,UOMTable$UOM_ID)]
dataset_uom <- conv_uom(dataset) 
dataset_uom<-dataset_uom%>%filter(!PEDIGREE_NAME=="FILLER")
dataset_uom<-dataset_uom%>%filter(!GROWSEASON%in%c("2015:08","2016:08","2017:08"))

################
# Diallel calculation Cases
# PCM3 x PCM3
# Hunch prediction
# Fall with summer testers
### files needed for next calculation
# gca_phenos_fall_PCM1_2021.Rdata subset of results only inbreds current cycle
# gca_tester_fall_PCM1_2021.Rdata subset of results only testers current cycle
# sumtst18_21_falleff_3yrTrain.Rdata  blup calculation using marker effects summer

##################
## PCM3 X PCM3

# Case 1 PCM3xPCM3
par_PCM3<-dataset_uom%>%filter(EXPER_STAGE_REF_ID=="P3")
par_PCM3_f<-na.omit(c(unique(unique(par_PCM3$P2)),unique(par_PCM3$P1)))

# select pedigrees complete diallel same test and inbred 
test<-unique(par_PCM3_f)
inbred<-unique(par_PCM3_f)

library(parallel)
#library(pbmcapply)
res<-mclapply(trait_choosen, function(x) diallel(x,results = res2,tester=test,inbr=inbred), mc.cores = 2,mc.preschedule = FALSE)
diallel_PCM3_Fall<-do.call(cbind,res)
diallel_PCM3_Fall <- diallel_PCM3_Fall[!(diallel_PCM3_Fall$P1 == diallel_PCM3_Fall$P2 ),]
#save(diallel_PCM3_Fall,file="/mnt/DH_processing/gca_fall 2021/diallel_PCM3_Fall.Rdata") ### now in AWS veg-apd-sdi-predictiveanalytcs-prod-workspace/ycao1/ITG/Cucumber/diallel_fall2021/diallel_PCM3_Fall.Rdata


###################################
## Case 2 hunch prediction inbreds with other testers from the same year
## see preparation file to create gca_phenos_fall_PCM1_2021.Rdata and gca_tester_fall_PCM1_2021.Rdata
# load GCA tables
#load("/mnt/DH_processing/gca_fall 2021/gca_phenos_fall_PCM1_2021.Rdata")
#load("/mnt/DH_processing/gca_fall 2021/gca_tester_fall_PCM1_2021.Rdata")
#s3save(gca_phenos_fall_PCM1_2021, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/gca_phenos_fall_PCM1_2021.Rdata")#training set all seasons
#s3save(gca_tester_fall_PCM1_2021, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/gca_tester_fall_PCM1_2021.Rdata")#training set all seasons

s3load(object = '/ycao1/ITG/Cucumber/diallel_fall2021/gca_phenos_fall_PCM1_2021.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

s3load(object = '/ycao1/ITG/Cucumber/diallel_fall2021/gca_tester_fall_PCM1_2021.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


tester = gca_tester_fall_PCM1_2021
inbred = gca_phenos_fall_PCM1_2021
res1<-rbind(inbred,tester)
# modify 
test<-unique(tester$PEDIGREE_NAME)
inbred<-unique(inbred$PEDIGREE_NAME)

res<-mclapply(trait, function(x) diallel(x,results = res1,tester=test,inbr=inbred), mc.cores = 2,mc.preschedule = FALSE)
diallel_DH2021_tst18_21<-do.call(cbind,res)
#save(diallel_DH2021_tst18_21,file="/mnt/DH_processing/gca_fall 2021/diallel_3yrDH2021_tst21.Rdata") ### INPUT FILE

###################################
## inbreds with summer testers 
#load("/mnt/DH_processing/gca_genos_sum2fall_tester/sumtst18_21_falleff_3yrTrain.Rdata") # marker effects from summer training set
#load("/mnt/DH_processing/gca_fall 2021/gca_phenos_fall_PCM1_2021.Rdata") # see preparation file
#s3save(df_fall, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/sumtst18_21_falleff_3yrTrain.Rdata")#training set all seasons

s3load(object = '/ycao1/ITG/Cucumber/diallel_fall2021/gca_phenos_fall_PCM1_2021.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
s3load(object = '/ycao1/ITG/Cucumber/diallel_fall2021/sumtst18_21_falleff_3yrTrain.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


# load GCA tables
tester = df_fall
inbred = gca_phenos_fall_PCM1_2021
inbred = inbred[,which(colnames(inbred)%in%colnames(tester))]
res1<-rbind(inbred,tester)

# modify 
test<-unique(tester$PEDIGREE_NAME)
inbred<-unique(inbred$PEDIGREE_NAME)

res<-mclapply(trait, function(x) diallel(x,results = res1,tester=test,inbr=inbred), mc.cores = 2,mc.preschedule = FALSE)
diallel_fallDH2021_sumtst18_21<-do.call(cbind,res)
#save(diallel_fallDH2021_sumtst18_21,file="/mnt/DH_processing/gca_fall 2021/diallel_fallDH2021_sumtst18_21.Rdata") ### INPUT FILE

