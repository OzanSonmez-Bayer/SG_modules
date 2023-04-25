library(readxl)
library(tidyverse)
library(aws.s3)
#load("/mnt/veg_blup/train_fall_2021/SSGBLUP_deep/AFW_C_deep/results_asreml.Rdata")#results file output reference for pedigrees

### results 
#load("/mnt/Cln_cucumber_fall09102021_new.Rdata")
s3load(object = '/ycao1/ITG/Cucumber/Fall/Cln_cucumber_fall12102021_new.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

trait_choosen<-c("AFW_C","EXTCO","FRLGT","FRNMK","SHAPE","WFSA_C","NFSA_C","FRNMKAvg_1-57","WFSA_CAvg_1-57","NFSA_CAvg_1-57","NFSA_CSum_1-57", "WFSA_CSum_1-57", "NETWTSum_1-57","LNUNI" )
path="/ycao1/ITG/Cucumber/GWS_results/veg_blup/train_fall_2021/SSGBLUP_deep/"
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

#####################################################
testers<-read_excel("/repos/ELBFA_GS_modules/geno_mod/Files_PCM/List of testers of PCM1 by season by year.xlsx",sheet = 1)
#check if all fall testers are calculated in results file
tst_fall_pcm1<-testers%>%filter(Season%in%c("fall") & Year %in% c("2020","2018","2019","2021"))
tst_fall_pcm1$Tester<-gsub("\\(E","Y3_E",tst_fall_pcm1$Tester)
tst_fall_pcm1$Tester<-gsub("\\)","",tst_fall_pcm1$Tester)
#table(unique(tst_fall_pcm1$Tester)%in%results$PEDIGREE_NAME)

#####################################################
###  tester pedigree gca from prediction SSGBLUP
#####################################################
dat0=data.frame(id=unique(tst_fall_pcm1$Tester)[which(unique(tst_fall_pcm1$Tester)%in%unique(results$PEDIGREE_NAME))])                           

dat=data.frame(id=unique(tst_fall_pcm1$Tester)[which(!unique(tst_fall_pcm1$Tester)%in%unique(results$PEDIGREE_NAME))])                           
dat$id<-gsub("((L","Y3_(L",dat$id,fixed = TRUE)
dat$id<-gsub("(LAUREEN:0111%","(LAUREEN:0111%)",dat$id,fixed = TRUE)
dat1<-data.frame(id=unique(dat$id)[which(unique(dat$id)%in%unique(results$PEDIGREE_NAME))])#      

#### missing
dat2<-data.frame(id=unique(dat$id)[which(!unique(dat$id)%in%unique(results$PEDIGREE_NAME))])                           
dat2$id<- gsub('.{4}$',"",dat2$id)
dat3<-data.frame(id=unique(dat2$id)[which(unique(dat2$id)%in%unique(results$PEDIGREE_NAME))])#

#missing replace
dat4<-data.frame(id=unique(dat2$id)[which(!unique(dat2$id)%in%unique(results$PEDIGREE_NAME))]) 
dat4$id<- gsub("Y3_EUR-Y316-7020GY/EUR-Y315-3038GY:0066%0001.0005.","Y3_EUR-Y316-7020GY/EUR-Y315-3038GY:0066%0001.",dat4$id)#
dat5<-data.frame(id=unique(dat4$id)[which(unique(dat4$id)%in%unique(results$PEDIGREE_NAME))]) 
table(dat4$id%in%results$PEDIGREE_NAME)
tst_fall_4sum<-rbind(dat0,dat1,dat3,dat4)

# match results with season data from pheno file
tst = pheno%>%filter(GROWSEASON=="2021:08")%>%distinct(PEDIGREE_NAME,.keep_all=TRUE)
tst = tst%>%filter(!PEDIGREE_NAME=="FILLER")
gca_tester_fall_PCM1_2021<-res2%>%filter(PEDIGREE_NAME%in%tst$P2)
table(unique(tst$P2)%in%res2$PEDIGREE_NAME)# 
unique(tst$P2)[unique(tst$P2)%in%res2$PEDIGREE_NAME]

#save(gca_tester_fall_PCM1_2021,file="/mnt/DH_processing/gca_fall 2021/gca_tester_fall_PCM1_2021.Rdata")### INPUT FILE
s3save(gca_tester_fall_PCM1_2021, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/gca_tester_fall_PCM1_2021.Rdata")#training set all seasons

############################
#### genos testers fall ####
############################

##### load inbred data
s3load(object = '/ycao1/ITG/Cucumber/inbreds/inbred_geno_cucumber.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

## training set 2020 GROWSEASON=="2021:08"
# DH lines processing_2020_fall
s3load(object = '/ycao1/ITG/Cucumber/inbreds/fall_DH_geno_cucumber.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

## training set 2021 GROWSEASON=="2020:08"
# DH lines processing_2021_fall
s3load(object = '/ycao1/ITG/Cucumber/inbreds/DH_2020_geno_cucumber.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#load("/mnt/DH_processing/DH_2020_geno_cucumber.Rdata")
#load("/mnt/DH_processing/fall_DH_geno_cucumber.Rdata")
#load("/mnt/inbred_geno_cucumber.Rdata")

names(DH_2020)[names(DH_2020) == 'ProgenyGermID'] <- 'GermID'
names(fall_DH)[names(fall_DH) == 'ProgenyGermID'] <- 'GermID'
all_inb<-rbind(inbred_geno,DH_2020,fall_DH)
dat_1<-data.frame(id=unique(tst_fall_4sum$id)[which(unique(tst_fall_4sum$id)%in%unique(all_inb$Pedigree))])#      

# missing no geno
table(tst_fall_4sum$id%in%all_inb$Pedigree)
dat_2<-data.frame(id=unique(tst_fall_4sum$id)[which(!unique(tst_fall_4sum$id)%in%unique(all_inb$Pedigree))])  

#Y3_E23L-2363:0178%0001.
#Y3_EUR-200-TEKST-GY/EURY314-3027GY:0038%0001.0001./EUR-200-BANC-GY:0060%0001.
#Y3_EUR-Y317-3057GY/EUR-200-BANC-GY:0056%0001.
#Y3_EUR-Y317-7027GY/EUR-Y317-5032GY:0049%0001.

genoDF_fall_4sum<-all_inb%>%filter(Pedigree%in%dat_1$id)%>%distinct(Pedigree,.keep_all=TRUE)
#save(genoDF_fall_4sum,file="/mnt/DH_processing/PCM1_DHFall_2021_genos/geno_testers18_20_fall.Rdata")### INPUT FILE
s3save(genoDF_fall_4sum, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/genoDF_fall_4sum.Rdata")#training set all seasons


############################
#### genos DH 2021 fall ####
############################

#extract phenos 2021 and subset P1 which are the DH Fall 2021

tst = pheno%>%filter(GROWSEASON=="2021:08")%>%distinct(PEDIGREE_NAME,.keep_all=TRUE)
tst = tst%>%filter(!PEDIGREE_NAME=="FILLER")
unique(tst$P1)

table(unique(tst$P1)%in%results$PEDIGREE_NAME)# 400 of 401 have gca
table(unique(tst$P1)%in%all_inb$Pedigree)# 317 - 84
miss = data.frame(id=unique(tst$P1)[!unique(tst$P1)%in%all_inb$Pedigree])
miss = drop_na(miss)

#Y3_EUR-Y318-1061GY/EUR-200-BANC-GY:0037%  Y3_EUR-Y318-1061GY/EUR-200-BANC-GY:0037%0001.
#Y3_EUR-Y318-1061GY/EUR-200-BANC-GY:0026%  Y3_EUR-Y318-1061GY/EUR-200-BANC-GY:0026%0001.
miss$id1<- gsub("JG1_","Y3_",miss$id)
miss$id1 = paste0(miss$id1,"0001.")
DH21_fall = miss%>%filter(id1%in%all_inb$Pedigree) 

table(unique(miss$id1)%in%all_inb$Pedigree)# 317 - 84
miss1 = miss%>%filter(!id1%in%all_inb$Pedigree) # no genos

DH_fall_PCM1_2021<-all_inb%>%filter(Pedigree%in%DH21_fall$id1)%>%distinct(Pedigree,.keep_all=TRUE)
#save(DH_fall_PCM1_2021,file="/mnt/DH_processing/PCM1_DHFall_2021_genos/DH_fall_PCM1_2021.Rdata")### INPUT FILE
s3save(DH_fall_PCM1_2021, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/DH_fall_PCM1_2021.Rdata")#training set all seasons


############################
#### gca DH 2021 fall ####
############################

tst = pheno%>%filter(GROWSEASON=="2021:08")%>%distinct(PEDIGREE_NAME,.keep_all=TRUE)
tst = tst%>%filter(!PEDIGREE_NAME=="FILLER")
gca_phenos_fall_PCM1_2021<-res2%>%filter(PEDIGREE_NAME%in%tst$P1)
table(unique(tst$P1)%in%res2$PEDIGREE_NAME)#1 missing
#save(gca_phenos_fall_PCM1_2021,file="/mnt/DH_processing/gca_fall 2021/gca_phenos_fall_PCM1_2021.Rdata") ### INPUT FILE
s3save(gca_phenos_fall_PCM1_2021, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/gca_phenos_fall_PCM1_2021.Rdata")#training set all seasons


