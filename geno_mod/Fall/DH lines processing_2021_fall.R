library(tidyverse)
library(readxl)

library(aws.s3)
library(plyr)
library(magrittr)
library(jsonlite)
library(data.table)
library(bit64)
library(tidyr)
library(httr)
library(readxl)
options(scipen = 999)


# Credentials

source('vaultCredentials.R')
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

geno_DH <- s3read_using(fread, 
                        object = "Cucumber/LongDutch/Y3/DH0_Spring_complete_20210914/HDgenotypes.txt.gz", 
                        bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')


length(unique(geno_DH$ProgenyGermID))

meta <- s3read_using(fread, 
                     object = "Cucumber/LongDutch/Y3/DH0_Spring_complete_20210914/Y3_LongDutch_2021SpringCycle_Combined_ProgenyTable.txt",
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

markerdata <- s3read_using(fread, 
                           object = 'Cucumber/LongDutch/Y3/DH0_Spring_complete_20210914/MarkerMetadata.txt',
                           bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

head(geno_DH)
head(meta)
head(markerdata)

# Wide formatting of the original genotypic file
geno_wide <- geno_DH %>% 
  dplyr::select(ProgenyGermID, MRN, HDgenotype) %>% 
  dcast(ProgenyGermID ~ MRN, value.var="HDgenotype")


FUN = function(x){
  unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
}


# Turn MRN to 0, 1, 2
#geno_wide_012 <- apply(geno_wide[, 2:ncol(geno_wide)], MARGIN = 2, 
#                       FUN = function(x){
#                         unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
#

#})

install.packages("future.apply")
library(future.apply)
library(tidyverse)
plan(multisession, workers = 10) 
geno_wide_012<-future_apply(geno_wide[, 2:ncol(geno_wide)],MARGIN = 2, FUN=FUN)

geno_wide<-data.frame(geno_wide)
geno_wide_012<-data.frame(geno_wide_012)
geno_wide_012 <- cbind(data.frame("ProgenyGermID"=geno_wide$ProgenyGermID), geno_wide_012[,1:ncol(geno_wide_012)])


# Pull Reference ID from S3

Crop <- 'Cucumber'
library(aws.s3)
obj <-get_object("s3://veg-apd-sdi-predictiveanalytcs-prod-reference-data/Cucumber/Cucumber_IDs.csv")  
csvcharobj <- rawToChar(obj)  
con <- textConnection(csvcharobj)  
CropIDs <- read.csv(file = con)

CropIDs_sub <- CropIDs %>%
  dplyr::select(M.GERMPLASM.X_ID,  M.GERMPLASM.PEDIGREE) %>%
  dplyr::rename(ProgenyGermID = M.GERMPLASM.X_ID,
                Pedigree = M.GERMPLASM.PEDIGREE) %>%
  mutate(ProgenyGermID = as.numeric(ProgenyGermID)) %>%
  unique()


geno_wide_012$ProgenyGermID <- as.numeric(geno_wide_012$ProgenyGermID)

geno_wide_ped <- CropIDs_sub %>% 
  dplyr::inner_join(geno_wide_012, by = 'ProgenyGermID')

dim(geno_wide_ped)
fall_DH<-geno_wide_ped
save(fall_DH,file="fall_DH_geno_cucumber.Rdata")


###########################################################

load("/mnt/DH_processing/fall_DH_geno_cucumber.Rdata")
load("/mnt/DH_processing/DH_2020_geno_cucumber.Rdata")

inbreds<-rbind(DH_2020,fall_DH)

### breeder reference
#dhSP<-read_excel("/mnt/DH_processing/List of testers of PCM1 by season by year.xlsx",sheet = 2)
dhSP<- s3read_using(FUN = read_excel, object = 's3://veg-apd-sdi-predictiveanalytcs-prod-workspace/ycao1/ITG/Cucumber/Files_PCM/Files_PCM/List of testers of PCM1 by season by year.xlsx', sheet = 2)

dhSP<-dhSP%>%dplyr::filter(Season=="fall")
dhSP<-dhSP%>%dplyr::filter(Year=="2021")# 
length(unique(dhSP$pedigree))
table(dhSP$Set_code)
#P1  P2 
#604  56 
table(dhSP$pedigree%in%CropIDs_sub$Pedigree)

dhSP<-dhSP%>%dplyr::filter(Set_code=="P1") #604 PCM crosses
dhSP<-dhSP[!duplicated(dhSP$Female_pedigree),] #393 Female pedigrees

#####
dhSP_1<-dhSP[which(dhSP$Female_pedigree%in%inbreds$Pedigree),]###########################
dhSP_1$Female_pedigree_1<-dhSP_1$Female_pedigree ########################### 4

table(dhSP_1$Female_pedigree%in%inbreds$Pedigree)#2

dhSP_2<-dhSP[which(!dhSP$Female_pedigree%in%dhSP_1$Female_pedigree),]
dhSP_2$Female_pedigree_1<- gsub("Y3_\\(","JG1_",dhSP_2$Female_pedigree)
dhSP_2$Female_pedigree_1<- gsub("0001)","",dhSP_2$Female_pedigree_1)
table(dhSP_2$Female_pedigree_1%in%inbreds$Pedigree)
dhSP_2<-dhSP_2[which(dhSP_2$Female_pedigree_1%in%inbreds$Pedigree),]####################### 130

#####
var<-c(unique(dhSP_2$Female_pedigree),unique(dhSP_1$Female_pedigree))
dhSP_3<-dhSP%>%dplyr::filter(!Female_pedigree%in%var)
dhSP_3<-dhSP_3[which(!dhSP_3$Female_pedigree%in%inbreds$Pedigree),]
dhSP_3$Female_pedigree_1<- gsub("Y3_\\(","Y3_",dhSP_3$Female_pedigree)
dhSP_3$Female_pedigree_1<- gsub("\\)","",dhSP_3$Female_pedigree_1)
dhSP_3$Female_pedigree_1<- gsub("%0001$","%0001.",dhSP_3$Female_pedigree_1)
dhSP_5<-dhSP_3[which(dhSP_3$Female_pedigree_1%in%inbreds$Pedigree),]############################### 205

#####
dhSP_4<-dhSP_3[which(!dhSP_3$Female_pedigree_1%in%inbreds$Pedigree),]
dhSP_4$Female_pedigree_1<- gsub("Y3_\\(","Y3_",dhSP_4$Female_pedigree)
dhSP_4$Female_pedigree_1<- gsub("%0001)","%0001.",dhSP_4$Female_pedigree_1)
inbreds$Pedigree[which(inbreds$Pedigree%in%dhSP_4$Female_pedigree_1)]
dhSP_6<-dhSP_4[which(dhSP_4$Female_pedigree_1%in%inbreds$Pedigree),]############################### 8

##### Not Found
dhSP_7<-dhSP_4[which(!dhSP_4$Female_pedigree_1%in%inbreds$Pedigree),]############################ 46

##### Using query file to match

ids_to_translate<-data.frame("peds"=inbreds$Pedigree)
save(ids_to_translate,file="DH_fall_genoids.Rdata")

dhSP$Female_pedigree_1<- gsub("Y3_\\(","Y3_",dhSP$Female_pedigree)
dhSP$Female_pedigree_1<- gsub("%0001)","%0001.",dhSP$Female_pedigree_1)

table(dhSP$Female_pedigree_1%in%inbredAll$OUTPUT_PEDIGREE_NAME)
dhSP_9<-dhSP[which(dhSP$Female_pedigree_1%in%inbredAll$OUTPUT_PEDIGREE_NAME),]



# not found query
dhSP_8<-dhSP[which(!dhSP$Female_pedigree_1%in%inbredAll$OUTPUT_PEDIGREE_NAME),]

######
## Found DH lines
all_fall<-rbind(dhSP_1,dhSP_2,dhSP_5,dhSP_6)
table(all_fall$Female_pedigree_1%in%inbreds$Pedigree)

######
## testers
######

table(unique(all_fall$`Male tester_pedigree`)%in%inbreds$Pedigree)
all_fall_1= all_fall[all_fall$`Male tester_pedigree`%in%inbreds$Pedigree,]
unique(as.character(all_fall$`Male tester_pedigree`))[which(!unique(all_fall$`Male tester_pedigree`)%in%inbreds$Pedigree)] 

#remove "(" and ) and last 4 digits
all_fall_2= all_fall[!all_fall$`Male tester_pedigree`%in%inbreds$Pedigree,]
all_fall_2$`Male tester_pedigree`<- gsub("\\(","Y3_",all_fall_2$`Male tester_pedigree`)
all_fall_2$`Male tester_pedigree`<- gsub('.{5}$',"",all_fall_2$`Male tester_pedigree`)
table(unique(all_fall_2$`Male tester_pedigree`)%in%inbreds$Pedigree)
all_fall_4= all_fall_2[!all_fall_2$`Male tester_pedigree`%in%inbreds$Pedigree,]

## keep only the one ones that match, one tester not found E23L-2363:0178%0001
all_fall_3= all_fall_2[all_fall_2$`Male tester_pedigree`%in%inbreds$Pedigree,]

dim(all_fall)

# final dataset after tester check
all_fall_final = rbind(all_fall_1,all_fall_3)#unique hybrids 295 out of 347 combinations missing due to tester E23L-2363:0178%0001

##### processing creating PCM1

fall_DH_germID<-inbreds[,1:2]

Dh_fall_1<-all_fall_final%>%left_join(fall_DH_germID,by=c("Female_pedigree_1"="Pedigree"))
dim(Dh_fall_1)
Dh_fall_1<-Dh_fall_1%>%select(Female_pedigree_1,`Male tester_pedigree`,ProgenyGermID,pedigree)%>%mutate(fem_id=ProgenyGermID)
Dh_fall_1<-Dh_fall_1%>%left_join(fall_DH_germID,by=c("Male tester_pedigree"="Pedigree"))
Dh_fall_1<-Dh_fall_1%>%select(Female_pedigree_1,`Male tester_pedigree`,ProgenyGermID.y,fem_id,pedigree)%>%mutate(test_id=ProgenyGermID.y)
Dh_fall_1<-Dh_fall_1%>%select(-ProgenyGermID.y)
head(Dh_fall_1)
# match Female_pedigree and Male tester pedigree
#@all the testcrosses in this subset match as expected
#for (i in 1:nrow(Dh_fall)) {
#  
#  Dh_fall[i,"Parent1GTStatus"] <- if_else(Dh_fall$Female_pedigree[i] %in% fall_DH$Pedigree,1,0)
#  Dh_fall[i,"Parent2GTStatus"] <- if_else(Dh_fall$`Male tester_pedigree`[i] %in% testSP_2021$Tester,1,0)
#  Dh_fall[i,"TotalParentsGT"] <- Dh_fall[i,"Parent1GTStatus"] + Dh_fall[i,"Parent2GTStatus"]
#}

#inbreds<-rbind(DH_2020,fall_DH)
#table(unique(Dh_fall_1$Female_pedigree_1)%in%inbreds$Pedigree)

hybridsWithBothParentsGT <- unique(Dh_fall_1$pedigree)#unique 295 crosses

### combine FALL and 2020 data
geno_DH <- s3read_using(fread, 
                        object = "Cucumber/LongDutch/Y3/DH0_Spring_complete_20210914/HDgenotypes.txt.gz", 
                        bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')



geno_DH_1 <- s3read_using(fread, 
                        object = "Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/HDgenotypes.txt.gz", 
                        bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
geno_DH=data.frame(geno_DH)
geno_DH_1=data.frame(geno_DH_1)

geno_DH=geno_DH[!as.numeric(geno_DH$ProgenyGermID)%in%as.numeric(geno_DH_1$ProgenyGermID),]
#geno_DH_1=geno_DH_1[!as.numeric(geno_DH_1$ProgenyGermID)%in%as.numeric(geno_DH$ProgenyGermID),]
table(as.numeric(geno_DH_1$ProgenyGermID)%in%as.numeric(geno_DH$ProgenyGermID))

all_geno_DH =rbind(geno_DH,geno_DH_1)

markerdata <- s3read_using(fread, 
                           object = 'Cucumber/LongDutch/Y3/DH0_Spring_complete_20210914/MarkerMetadata.txt',
                           bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')


####

finalFileName <- "/mnt/Y3_fallSeasonPCMSet2021_SynthHybrids.tab"

for (h in 1:nrow(Dh_fall_1)) {
  print(paste0("Currently working on hybrid number ",h," of ", nrow(Dh_fall_1)))

  HybridTable <- data.frame(HybridGermID=rep(hybridsWithBothParentsGT[h],nrow(markerdata)),MRN=markerdata$MRN) %>% distinct()      
  
  bothParentsImputeGT <- data.frame(MRN=markerdata$MRN)
  
  currentParent1 <- Dh_fall_1[h,"fem_id"]
  
  all_geno_DH<-data.frame(all_geno_DH)
  all_geno_DH$ProgenyGermID<-as.numeric(all_geno_DH$ProgenyGermID)
  
  fem<-all_geno_DH[all_geno_DH$ProgenyGermID%in%currentParent1$fem_id,]
  
  currentParent2 <- Dh_fall_1[h,"test_id"]
  tester<-all_geno_DH[all_geno_DH$ProgenyGermID%in%currentParent2$test_id,]
  
  bothParentsImputeGT <- merge(bothParentsImputeGT,fem[,c(2,5)],all.x = T)
  colnames(bothParentsImputeGT) <- c("MRN","Parent1")
  bothParentsImputeGT <- merge(bothParentsImputeGT,tester[c(2,5)], all.x = T)
  colnames(bothParentsImputeGT) <- c("MRN","Parent1","Parent2")
  
  
  synthesizedHybridTable <- data.frame(HybridGermID=NA,MRN=NA,Genotype=NA,Probability=NA) %>% distinct()                               
  
  parentCondition1 <- bothParentsImputeGT[which(bothParentsImputeGT$Parent1==bothParentsImputeGT$Parent2),]
  
  synthHybrid_condition1 <- merge(HybridTable,parentCondition1[,1:2],by = "MRN")                                         
  
  
  if (nrow(synthHybrid_condition1) >= 1) {
    synthHybrid_condition1$Probability <- 1.0
    colnames(synthHybrid_condition1) <- c("MRN","HybridGermID","Genotype","Probability")
    synthHybrid_condition1 <- synthHybrid_condition1[,c("HybridGermID","MRN","Genotype","Probability")]
    synthesizedHybridTable <- rbind(synthesizedHybridTable,synthHybrid_condition1)
  }
  
  ## Condition 2: Parent1 is 00 and Parent 2 is 11
  parentCondition2 <- bothParentsImputeGT[which(bothParentsImputeGT$Parent1 == "0|0" & bothParentsImputeGT$Parent2 == "1|1"),]
  synthHybrid_condition2 <- merge(HybridTable,parentCondition2[,1:2],by = "MRN")
  if (nrow(synthHybrid_condition2) >= 1) {  
    synthHybrid_condition2$Probability <- 1.0
    colnames(synthHybrid_condition2) <- c("MRN","HybridGermID","Genotype","Probability")
    synthHybrid_condition2$Genotype <- "0|1"
    synthHybrid_condition2 <- synthHybrid_condition2[,c("HybridGermID","MRN","Genotype","Probability")]
    synthesizedHybridTable <- rbind(synthesizedHybridTable,synthHybrid_condition2)
  }
  
  ## Condition 3: Parent1 is 11 and Parent 2 is 00
  parentCondition3 <- bothParentsImputeGT[which(bothParentsImputeGT$Parent1 == "1|1" & bothParentsImputeGT$Parent2 == "0|0"),]
  synthHybrid_condition3 <- merge(HybridTable,parentCondition3[,1:2],by = "MRN")
  if(nrow(synthHybrid_condition3) >= 1){
    synthHybrid_condition3$Probability <- 1.0
    colnames(synthHybrid_condition3) <- c("MRN","HybridGermID","Genotype","Probability")
    synthHybrid_condition3$Genotype <- "1|0"
    synthHybrid_condition3 <- synthHybrid_condition3[,c("HybridGermID","MRN","Genotype","Probability")]
    synthesizedHybridTable <- rbind(synthesizedHybridTable,synthHybrid_condition3)
  }
  
  ## Condition 4: Parent 1 is 01 and Parent 2 is 10
  parentCondition4 <- bothParentsImputeGT[which(bothParentsImputeGT$Parent1 == "0|1" & bothParentsImputeGT$Parent2 == "1|0"),]
  synthHybrid_condition4 <- merge(HybridTable,parentCondition4[,1:2],by = "MRN")
  if (nrow(synthHybrid_condition4) >= 1) {
    synthHybrid_condition4$Probability <- 0.5
    colnames(synthHybrid_condition4) <- c("MRN","HybridGermID","Genotype","Probability")
    synthHybrid_condition4$Genotype <- "1|0"
    synthHybrid_condition4 <- synthHybrid_condition4[,c("HybridGermID","MRN","Genotype","Probability")]
    synthesizedHybridTable <- rbind(synthesizedHybridTable,synthHybrid_condition4)
  }
  
  ## Condition 5: Parent 1 is 01 and Parent 2 is 01 
  parentCondition5 <- bothParentsImputeGT[which(bothParentsImputeGT$Parent1 == "1|0" & bothParentsImputeGT$Parent2 == "0|1"),]
  synthHybrid_condition5 <- merge(HybridTable,parentCondition5[,1:2],by = "MRN")
  if (nrow(synthHybrid_condition5) >= 1) {
    synthHybrid_condition5$Probability <- 0.5
    colnames(synthHybrid_condition5) <- c("MRN","HybridGermID","Genotype","Probability")
    synthHybrid_condition5$Genotype <- "1|0"
    synthHybrid_condition5 <- synthHybrid_condition5[,c("HybridGermID","MRN","Genotype","Probability")]
    synthesizedHybridTable <- rbind(synthesizedHybridTable,synthHybrid_condition5)
  }
  
  ## Condition 6: 1 Parent Het and 1 Parent Hom
  parentCondition6 <- bothParentsImputeGT[!(bothParentsImputeGT$MRN %in% synthesizedHybridTable$MRN),]
  synthHybrid_condition6 <- merge(HybridTable,parentCondition6[,1:2],by = "MRN")
  if (nrow(synthHybrid_condition6) >= 1) {
    synthHybrid_condition6$Probability <- 0.5
    colnames(synthHybrid_condition6) <- c("MRN","HybridGermID","Genotype","Probability")
    synthHybrid_condition6$Genotype <- "NA"
    synthHybrid_condition6 <- synthHybrid_condition6[,c("HybridGermID","MRN","Genotype","Probability")]
    synthesizedHybridTable <- rbind(synthesizedHybridTable,synthHybrid_condition6)
  }
  
  synthesizedHybridTable<-synthesizedHybridTable%>%drop_na(HybridGermID)
  
  write.table(synthesizedHybridTable,file = finalFileName,append = file.exists(finalFileName),sep = "\t",row.names = F,col.names=!file.exists(finalFileName)) 
  #write.table(synthesizedHybridTable,file = paste0(getwd(),"/IndividualSynthHybridFiles/",currentHybrid,"_SynthGeno.tab"),append = FALSE,sep = "\t",row.names = F)
}



# Wide formatting of the original genotypic file
#hybrid_geno <- read.csv(file = "/mnt/Y3_fallSeasonPCMSet2021_SynthHybrids_1.tab", sep="\t", header=T)
hybrid_geno <- read.csv(file = "/mnt/Y3_fallSeasonPCMSet2021_SynthHybrids.tab", sep="\t", header=T)

geno_wide <- hybrid_geno %>% 
  dplyr::select(HybridGermID, MRN, Genotype) %>% 
  dcast(HybridGermID ~ MRN, value.var="Genotype")


FUN = function(x){
  unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
}


# Turn MRN to 0, 1, 2
#geno_wide_012 <- apply(geno_wide[, 2:ncol(geno_wide)], MARGIN = 2, 
#                       FUN = function(x){
#                         unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
#

#})

install.packages("future.apply")
library(future.apply)
plan(multisession, workers = 10) 
geno_wide_012<-future_apply(geno_wide[, 2:ncol(geno_wide)],MARGIN = 2, FUN=FUN)

geno_wide<-data.frame(geno_wide)
geno_wide_012<-data.frame(geno_wide_012)
geno_wide_012 <- cbind(data.frame("ProgenyGermID"=geno_wide$HybridGermID), geno_wide_012[,1:ncol(geno_wide_012)])


# Pull Reference ID from S3

Crop <- 'Cucumber'
library(aws.s3)
obj <-get_object("s3://veg-apd-sdi-predictiveanalytcs-prod-reference-data/Cucumber/Cucumber_IDs.csv")  
csvcharobj <- rawToChar(obj)  
con <- textConnection(csvcharobj)  
CropIDs <- read.csv(file = con)

CropIDs_sub <- CropIDs %>%
  dplyr::select(M.GERMPLASM.X_ID,  M.GERMPLASM.PEDIGREE) %>%
  dplyr::rename(ProgenyGermID = M.GERMPLASM.X_ID,
                Pedigree = M.GERMPLASM.PEDIGREE) %>%
  mutate(ProgenyGermID = as.numeric(ProgenyGermID)) %>%
  unique()


geno_wide_ped <- CropIDs_sub %>% 
  dplyr::inner_join(geno_wide_012, by = c('Pedigree'='ProgenyGermID'))

dim(geno_wide_ped)

Fall_21_train<-geno_wide_ped
save(Fall_21_train,file="/mnt/DH_processing/Fall_2021_geno_traincuc.Rdata")
