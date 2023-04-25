

library(aws.s3)
library(plyr)
library(magrittr)
library(jsonlite)
library(data.table)
library(bit64)
library(tidyr)
library(httr)

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
                     object = "Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/HDgenotypes.txt.gz", 
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')


length(unique(geno_DH$ProgenyGermID))#2783

meta <- s3read_using(fread, 
                     object = "Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/WorkingFiles/Cucumber_LongDutch_Y3_2500Lines_ProgenyTable.txt",
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

markerdata <- s3read_using(fread, 
                           object = 'Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/MarkerMetadata.txt',
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
DH_2020<-geno_wide_ped
save(DH_2020,file="DH_2020_geno_cucumber.Rdata")
#save(summer_DH,file="summer_DH_geno_cucumber.Rdata")


###### all DH in PCM1 crosses file
library(readxl)
#dhSP_inv<-read_excel("/mnt/DH_processing/List of testers of PCM1 by season by year_InBID.xlsx",sheet = 2)
dhSP_inv<- s3read_using(FUN = read_excel, object = 's3://veg-apd-sdi-predictiveanalytcs-prod-workspace/ycao1/ITG/Cucumber/Files_PCM/Files_PCM/List of testers of PCM1 by season by year_InBID.xlsx', sheet = 2)

#dhSP<-read_excel("/mnt/List of testers of PCM1 by season by year.xlsx",sheet = 2)
dhSP<- s3read_using(FUN = read_excel, object = 's3://veg-apd-sdi-predictiveanalytcs-prod-workspace/ycao1/ITG/Cucumber/Files_PCM/Files_PCM/List of testers of PCM1 by season by year.xlsx', sheet = 2)

#dhSP<-dhSP%>%dplyr::filter(Season=="summer")
dhSP<-dhSP%>%dplyr::filter(Year=="2020")
dhSP_1<-dhSP%>%dplyr::filter(Set_code=="P1")
dhSP_2<-dhSP%>%dplyr::filter(Set_code=="P2")


######################################
###### hybrids expected genomics team

dat<-fread("/mnt/Y3_HybridInvPedTable.tab")
dat_2020<-dat%>%filter(HybridSeason%in%c("2020:05","2020:08"))#1708 hybrids
dat_2020_05<-dat%>%filter(HybridSeason%in%c("2020:05"))#1708 hybrids
dat_2020_08<-dat%>%filter(HybridSeason%in%c("2020:08"))#1708 hybrids

#############################
#create dat_2020 hybrids fall
#############################

load("/mnt/Cln_cucumber_fall06102021_new.Rdata")
phen_2020_fall<-pheno%>%filter(GROWSEASON=="2020:08")
table(unique(phen_2020_fall$PEDIGREE_NAME)%in%unique(dat_2020_08$HybridPedigree))#513 T 27 F
length(unique(dat_2020_08$HybridPedigree))

FL<-dat_2020_08%>%separate(HybridPedigree,into=c("A","B"),sep="\\+")
FL$HybridPedigree<-dat_2020_08$HybridPedigree
FL_1<-FL[FL$A%in%DH_2020$Pedigree,] #
FL_1$Female_pedigree<-FL_1$A #11

FL_2<-FL[!FL$A%in%DH_2020$Pedigree,]
FL_2$Female_pedigree<- gsub("Y3_\\(","Y3_",FL_2$A)
FL_2$Female_pedigree<- gsub("\\)","",FL_2$Female_pedigree)
FL_2$Female_pedigree<- gsub("%0001$","%0001.",FL_2$Female_pedigree)
FL_3<-FL_2[FL_2$Female_pedigree%in%DH_2020$Pedigree,] #346
table(FL_3$Female_pedigree%in%DH_2020$Pedigree)


FL_4<-FL_2[!FL_2$Female_pedigree%in%DH_2020$Pedigree,]
FL_4$Female_pedigree<- gsub("Y3_\\(","Y3_",FL_4$A)
FL_4$Female_pedigree<- gsub("%0001)","%0001.",FL_4$Female_pedigree)
FL_5<-FL_4[FL_4$Female_pedigree%in%DH_2020$Pedigree,]#22
table(FL_5$Female_pedigree%in%DH_2020$Pedigree)


FL_6<-FL_4[!FL_4$Female_pedigree%in%DH_2020$Pedigree,] 
FL_6$Female_pedigree<- gsub("\\)","",FL_6$Female_pedigree)#108 peds
FL_6$Female_pedigree<-gsub('.{4}$', '',FL_6$Female_pedigree)
table(unique(FL_6$Female_pedigree)%in%DH_2020$Pedigree)
FL_7<-FL_6[FL_6$Female_pedigree%in%DH_2020$Pedigree,]#41


FL_9<-FL_6[!FL_6$Female_pedigree%in%DH_2020$Pedigree,]
FL_9$Female_pedigree<- gsub("Y3_\\(","Y3_",FL_9$A)
FL_9$Female_pedigree<-gsub('.{1}$', '',FL_9$Female_pedigree)
FL_9$Female_pedigree<-gsub('.{4}$', '',FL_9$Female_pedigree)
FL_10<-FL_9[FL_9$Female_pedigree%in%DH_2020$Pedigree,] #20

table(unique(FL_10$Female_pedigree)%in%DH_2020$Pedigree)

FL_11<-FL_9[!FL_9$Female_pedigree%in%DH_2020$Pedigree,]
load("/mnt/sql_2020_train_genos.Rdata")
table(FL_11$A%in%inbredAll$OUTPUT_PEDIGREE_NAME)
FL_12 <- FL_11[FL_11$Female_pedigree%in%inbredAll$OUTPUT_PEDIGREE_NAME,]
FL_13 <-FL_12%>%left_join(inbredAll,by= c("Female_pedigree"="OUTPUT_PEDIGREE_NAME"))
FL_13$Female_pedigree<-FL_13$INPUT_PEDIGREE_NAME
FL_13<-FL_13[!duplicated(FL_13$Female_pedigree),]
FL_13<-data.frame(FL_13)
FL_13<-FL_13[,which(colnames(FL_13)%in%colnames(FL_11))] #9

####
FL_14 <- FL_11[!FL_11$Female_pedigree%in%inbredAll$OUTPUT_PEDIGREE_NAME,]#not found
length(unique(FL_14$Female_pedigree))#22
length(unique(FL_14$A))#28

### all data
all_fall_20<-rbind(FL_1,FL_3,FL_5,FL_7,FL_10,FL_13)
length(unique(all_fall_20$Female_pedigree))#449
length(unique(all_fall_20$A))#457

#tester
table(unique(all_fall_20$B)%in%DH_2020$Pedigree)#17, 7 no
load("/mnt/inbred_geno_cucumber.Rdata")
unique(all_fall_20$B)[which(!unique(all_fall_20$B)%in%DH_2020$Pedigree)]
table(unique(all_fall_20$B)%in%inbred_geno$Pedigree)#17, 7 not in inbred file previously delivered
unique(all_fall_20$B)[which(!unique(all_fall_20$B)%in%inbred_geno$Pedigree)]

FL_T1<-all_fall_20[all_fall_20$B%in%DH_2020$Pedigree,]
FL_T2<-all_fall_20[!all_fall_20$B%in%DH_2020$Pedigree,]
FL_T2$B<- gsub("\\(","",FL_T2$B)
FL_T2$B<- gsub("\\)","",FL_T2$B)
FL_T2$B<-gsub('.{4}$', '',FL_T2$B)
FL_T2$B<-paste0("Y3_",FL_T2$B)
table(unique(FL_T2$B)%in%DH_2020$Pedigree)
unique(FL_T2$B)[which(!unique(FL_T2$B)%in%DH_2020$Pedigree)]

#Y3_LAUREEN:0111%/EURY314-1018GY:0080%0001.0001.   Y3_(LAUREEN:0111%)/EURY314-1018GY:0080%0001.0001.
#Y3_EUR-200-TEKST-GY/EURY314-3027GY:0038%0001.0001./EUR-200-BANC-GY:0060%0001. to Y3_(EUR-200-TEKST-GY/EURY314-3027GY:0038%0001.0001.)/EUR-200-BANC-GY:0060%0001.	
FL_T2$B[FL_T2$B=="Y3_LAUREEN:0111%/EURY314-1018GY:0080%0001.0001."] <- "Y3_(LAUREEN:0111%)/EURY314-1018GY:0080%0001.0001."
FL_T2$B[FL_T2$B=="Y3_EUR-200-TEKST-GY/EURY314-3027GY:0038%0001.0001./EUR-200-BANC-GY:0060%0001."] <- "Y3_(EUR-200-TEKST-GY/EURY314-3027GY:0038%0001.0001.)/EUR-200-BANC-GY:0060%0001."

FL_alldata=rbind(FL_T1,FL_T2)
FL_alldata$Tester<-FL_alldata$B

### Get progenyID to create crosses
DH_germID_20<-DH_2020[,1:2]
Dh_fall_1<-FL_alldata%>%left_join(DH_germID_20,by=c("Female_pedigree"="Pedigree"))
Dh_fall_1<-Dh_fall_1%>%select(Female_pedigree,Tester,ProgenyGermID,HybridPedigree)%>%mutate(fem_id=ProgenyGermID)

Dh_fall_1<-Dh_fall_1%>%left_join(DH_germID_20,by=c("Tester"="Pedigree"))
colnames(Dh_fall_1)
Dh_fall_1<-Dh_fall_1%>%select(Female_pedigree,Tester,ProgenyGermID.y,fem_id,HybridPedigree)%>%mutate(test_id=ProgenyGermID.y)
Dh_fall_1<-Dh_fall_1%>%select(-ProgenyGermID.y)


#create dat_2020 hybrids summer
#load("/mnt/Cln_cucumber_summer27102021_new.Rdata")
#phen_2020_summer<-pheno_summer%>%filter(GROWSEASON=="2020:05")
#table(unique(phen_2020_summer$PEDIGREE_NAME)%in%unique(dat_2020_05$HybridPedigree))#513 T 27 F

hybridsWithBothParentsGT <- unique(Dh_fall_1$HybridPedigree)

####

finalFileName <- "/mnt/Y3_DH_fall_2020_SynthHybrids.tab"

for (h in 1:nrow(Dh_fall_1)) {

  print(paste0("Currently working on hybrid number ",h," of ", nrow(Dh_fall_1)))

  
HybridTable <- data.frame(HybridGermID=rep(hybridsWithBothParentsGT[h],nrow(markerdata)),MRN=markerdata$MRN) %>% distinct()      

bothParentsImputeGT <- data.frame(MRN=markerdata$MRN)

currentParent1 <- Dh_fall_1[h,"fem_id"]

geno_DH<-data.frame(geno_DH)
fem<-geno_DH[geno_DH$ProgenyGermID%in%currentParent1$fem_id,]

currentParent2 <- Dh_fall_1[h,"test_id"]
tester<-geno_DH[geno_DH$ProgenyGermID%in%currentParent2$test_id,]

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
synthesizedHybridTable<-synthesizedHybridTable%>%drop_na(Genotype)

write.table(synthesizedHybridTable,file = finalFileName,append = file.exists(finalFileName),sep = "\t",row.names = F,col.names=!file.exists(finalFileName)) 
#write.table(synthesizedHybridTable,file = paste0(getwd(),"/IndividualSynthHybridFiles/",currentHybrid,"_SynthGeno.tab"),append = FALSE,sep = "\t",row.names = F)
}


#####################
## Load spring PCM 1 Spring
#bothParentsImputeGT[which(!bothParentsImputeGT$MRN%in%synthesizedHybridTable$MRN),]
hybrid_geno <- read.csv(file = "/mnt/Y3_DH_fall_2020_SynthHybrids.tab", sep="\t", header=T)

hybrid_wide <- hybrid_geno %>% 
  dplyr::select(HybridGermID, MRN, Genotype) %>% 
  dcast(HybridGermID ~ MRN, value.var="Genotype")#491

FUN = function(x){
  unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
}

library(future)
plan(multisession, workers = 10) 
geno_wide_012<-future_apply(hybrid_wide[, 2:ncol(hybrid_wide)],MARGIN = 2, FUN=FUN)
geno_wide_012 <- cbind(data.frame("ProgenyGermID"=hybrid_wide$HybridGermID), geno_wide_012[,1:ncol(geno_wide_012)])

# Check to see if there are any NAs in hybrid geno file

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
Fall_20_train<-geno_wide_ped
save(Fall_20_train,file="Fall_2020_geno_traincuc.Rdata")
