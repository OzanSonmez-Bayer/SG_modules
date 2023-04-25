##############
## results were additionally formatted to add ASORT1
##############
# loaded to AWS diallel fall 2021
#/mnt/DH_processing/gca_fall 2021/diallel_3yrDH2021_tst21.Rdata
#/mnt/DH_processing/gca_fall 2021/diallel_fallDH2021_sumtst18_21.Rdata
#/mnt/DH_processing/gca_fall 2021/diallel_PCM3_Fall.Rdata
#"diallel_DH2021_tst18_21"        "diallel_fallDH2021_sumtst18_21" "diallel_PCM3_Fall"   
#s3save(diallel_DH2021_tst18_21, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/diallel_3yrDH2021_tst21.Rdata")#training set all seasons
#s3save(diallel_fallDH2021_sumtst18_21, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/diallel_fallDH2021_sumtst18_21.Rdata")#training set all seasons
#s3save(diallel_PCM3_Fall, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = "/ycao1/ITG/Cucumber/diallel_fall2021/diallel_PCM3_Fall.Rdata")#training set all seasons

library(readxl)
library(dplyr)
library(tidyverse)
library(magrittr)

################
### PCM1
test=read_xlsx("/repos/ELBFA_GS_modules/geno_mod/Files_PCM/List of testers of PCM1 by season by year_InBID_1.xlsx",sheet = 6)
s3load(object = "/ycao1/ITG/Cucumber/diallel_fall2021/diallel_3yrDH2021_tst21.Rdata", 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
#load("/mnt/DH_processing/gca_fall 2021/diallel_3yrDH2021_tst21.Rdata")
dat<-diallel_DH2021_tst18_21
ref=dat[,c("P1","P2")]
ref$P2_1=ref$P2
ref$P2_1=gsub("JG1_","Y3_",ref$P2_1)
ref1=ref[ref$P2_1%in%test$Pedigree,]
ref2=ref[!ref$P2_1%in%test$Pedigree,]
ref2$P2_1<-paste0(ref2$P2_1,"0001.")
all_ref=rbind(ref1,ref2)
table(unique(all_ref$P2_1)%in%test$Pedigree) #387

test=test[,c("Pedigree","ASRT1...1")]
df_to_join1 <- unique(test)

all_ref=left_join(all_ref,df_to_join1,by=c("P2_1"="Pedigree"))
all_ref=left_join(all_ref,df_to_join1,by=c("P1"="Pedigree"))

P1ref=all_ref[,c("P1","ASRT1...1.y")]
P2ref=all_ref[,c("P2","P2_1","ASRT1...1.x")]

P1ref=P1ref[!duplicated(P1ref$P1),]
P2ref=P2ref[!duplicated(P2ref$P2),]

dat_asort= left_join(P1ref,dat,by="P1")
dat_asort=left_join(P2ref,dat_asort,by="P2")
dat_asort$ASTR1_P1=dat_asort$ASRT1...1.y
dat_asort$ASTR1_P2=dat_asort$ASRT1...1.x
#write.table(dat_asort,file="C:/Users/ELBFA/Downloads/cucumber/diallel_3yrDH2021_tst21_asort1.csv",col.names = TRUE,row.names = FALSE,quote = FALSE,sep ="," )

#################
#### Fall with summer
test=read_xlsx("/repos/ELBFA_GS_modules/geno_mod/Files_PCM/List of testers of PCM1 by season by year_InBID_1.xlsx",sheet = 6)
#load(file="/mnt/DH_processing/gca_fall 2021/diallel_fallDH2021_sumtst18_21.Rdata")
s3load(object = "/ycao1/ITG/Cucumber/diallel_fall2021/diallel_fallDH2021_sumtst18_21.Rdata", 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
dat2<-diallel_fallDH2021_sumtst18_21
ref=dat2[,c("P1","P2")]
ref$P2_1=ref$P2
ref$P2_1=gsub("JG1_","Y3_",ref$P2_1)
ref1=ref[ref$P2_1%in%test$Pedigree,]
ref2=ref[!ref$P2_1%in%test$Pedigree,]
ref2$P2_1<-paste0(ref2$P2_1,"0001.")
all_ref=rbind(ref1,ref2)
table(unique(all_ref$P2_1)%in%test$Pedigree) #387
table(unique(all_ref$P1)%in%test$Pedigree) #387
test=test[,c("Pedigree","ASRT1...1")]
df_to_join1 <- unique(test)

all_ref=left_join(all_ref,df_to_join1,by=c("P2_1"="Pedigree"))
all_ref=left_join(all_ref,df_to_join1,by=c("P1"="Pedigree"))

P1ref=all_ref[,c("P1","ASRT1...1.y")]
P2ref=all_ref[,c("P2","P2_1","ASRT1...1.x")]

P1ref=P1ref[!duplicated(P1ref$P1),]
P2ref=P2ref[!duplicated(P2ref$P2),]

dat_asort= left_join(P1ref,dat2,by="P1")
dat_asort=left_join(P2ref,dat_asort,by="P2")
dat_asort$ASTR1_P1=dat_asort$ASRT1...1.y
dat_asort$ASTR1_P2=dat_asort$ASRT1...1.x
#write.table(dat_asort,file="diallel_fallDH2021_sumtst18_21_asort1.csv",col.names = TRUE,row.names = FALSE,quote = FALSE,sep ="," )

#################
#################
###### CxC
test=read_xlsx("/repos/ELBFA_GS_modules/geno_mod/Files_PCM/List of testers of PCM1 by season by year_InBID_1.xlsx",sheet = 6)
#load(file="/mnt/DH_processing/gca_fall 2021/diallel_PCM3_Fall.Rdata")
s3load(object = "/ycao1/ITG/Cucumber/diallel_fall2021/diallel_PCM3_Fall.Rdata", 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
dat3<-diallel_PCM3_Fall
subs=test[,colnames(test)%in%c("Pedigree","ASRT1...1")]
df_to_join <- unique(subs)
dat_pcm3=inner_join(dat3,df_to_join,by=c("P1"="Pedigree"))
dat_pcm3=inner_join(dat_pcm3,df_to_join,by=c("P2"="Pedigree"))
#write.table(dat_pcm3,file="fall_pcm3_cxc.csv",col.names = TRUE,row.names = FALSE,quote = FALSE,sep ="," )

