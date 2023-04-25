#### libraries
library(readxl)
library(tidyverse)

##################
# 2022 summer
##################

dhSP<-read_excel("/mnt/DH_processing/List of hybrids to be tested in Summer 2021-22.xlsx",sheet = "summer 2022")
miss_inv5=data.frame(id=unique(c(dhSP$`Inv BID_Female`,dhSP$`Inv BID_Male`)))# input sql function
miss_inv5=drop_na(miss_inv5)

#miss_inv5_1=data.frame(id=dhSP3$`Inv BID`)# input sql function
miss_inv3=data.frame(id=unique(c(dhSP$Pedigree_Female,dhSP$Pedigree_Male)))# input sql function
miss_inv3=drop_na(miss_inv3)

dhSP=data.frame(Inv_BID=c(dhSP$`Inv BID_Female`,dhSP$`Inv BID_Male`),Pedigree=c(dhSP$Pedigree_Female,dhSP$Pedigree_Male))

#### Parents Pairs
dhSP1<-read_excel("/mnt/DH_processing/List of hybrids to be tested in Summer 2021-22.xlsx",sheet = "summer 2022")
dhSP1<-dhSP1%>%select(Pedigree,`Inv BID_Female`,`Inv BID_Male`)

#sql
#x1=as.character(miss_inv3$id)
#max=10
#y1=seq_along(miss_inv3$id)
#invs1 <- split(x1, ceiling(y1/max))
invs1=list("Y3_EUR-Y316-3049GY/EUR-Y315-7014GY:0064%0001.")
Allinbred1=ped_name_sql(invs1)
#save(Allinbred1,file="/mnt/DH_processing/sql_pedigree_name_summer2020.Rdata")


#x=as.character(miss_inv5$id)
#max=10
#y=seq_along(miss_inv5$id)
#invs <- split(x, ceiling(y/max))
#Allinbred3=bar_inv_sql(invs)
#save(Allinbred3,file="/mnt/DH_processing/sql_invbid_summer2020.Rdata")

#### load genotypes
load("/mnt/DH_processing/DH_2020_geno_cucumber.Rdata")#2020
load("/mnt/inbred_geno_cucumber.Rdata")#training
load("/mnt/DH_processing/summer_DH_geno_cucumber.Rdata")#2021

#### Geno format
inbred_geno=dplyr::rename(inbred_geno,c("ProgenyGermID"="GermID") )
inb=rbind(summer_DH,DH_2020,inbred_geno)

load("/mnt/DH_processing/sql_pedigree_name_summer2022.Rdata")
load("/mnt/DH_processing/sql_invbid_summer2022.Rdata")

############################
#### input geno automat ####
############################
miss_inv5=miss_inv5
inbredAll=Allinbred3
inbredAll2=Allinbred1
inb=inb
dhSP1=dhSP1

ids_par<-function(miss_inv5,inbredAll,inbredAll2,inb,dhSP1){
  # @inbredAll sql file matched by Inv BID/Barcode
  # @inbredAll2 sql file matched by Pedigree Name
  # @inb genotyped data combined
  # @dhSP1: file =  "Pedigree" "InvBID_Female" "InvBID_Male"  
  # @miss_inv5 Inv BID of parents that was input to sql
  
  Dh_fall_2<-miss_inv5%>%left_join(inbredAll,by=c("id"="INPUT_BARCODE"))
  
  ## found input germplasm
  geno1=Dh_fall_2[Dh_fall_2$INPUT_GERMPLASM_ID%in%inb$ProgenyGermID,]#use 131 ID
  data1=data.frame(id=geno1$id,GermID=geno1$INPUT_GERMPLASM_ID)
  unique(data1$id)#397
  
  ## found output germplasm
  geno2=Dh_fall_2[!Dh_fall_2$INPUT_GERMPLASM_ID%in%inb$ProgenyGermID,]
  geno3=geno2[as.numeric(geno2$OUTPUT_GERMPLASM_ID)%in%as.numeric(inb$ProgenyGermID),]
  data2=data.frame(id=geno3$id,GermID=geno3$OUTPUT_GERMPLASM_ID)#use 17 ID
  table(unique(as.character(miss_inv5$id))%in%c(as.character(data1$id),as.character(data2$id)))#False
  
  if(dim(data2)[1]==0){
    alldata=data1
    
  }else{
    alldata=rbind(data1,data2)
  }
  
  ####################################################
  #create a file for matching geno file and phenofile
  ####################################################
  
  #missing 
  data3=miss_inv5%>%dplyr::filter(!id%in%alldata$id)
  data3=data3%>%left_join(dhSP,by=c("id"="Inv_BID"))%>%distinct()
  data3=data3%>%left_join(Allinbred1,by=c("Pedigree"="INPUT_PEDIGREE_NAME"))%>%drop_na()%>%select("id","Pedigree","INPUT_GERMPLASM_ID")%>%distinct()
  data3=data3%>%left_join(inb)%>%drop_na()
  data3=data.frame(id=data3$id,GermID=data3$ProgenyGermID)
  
  if(dim(data3)[1]==0){
    alldata=alldata
  }else{
    alldata=rbind(alldata,data3)
  }
  
  miss_seas=miss_inv5%>%dplyr::filter(!id%in%alldata$id)
  miss_seas=miss_seas%>%left_join(inbredAll,by=c("id"="INPUT_BARCODE"))%>%select(id,INPUT_GERMPLASM_ID,OUTPUT_GERMPLASM_ID,INPUT_PEDIGREE_NAME)%>%distinct()
  miss1=miss_seas%>%left_join(dhSP,by=c("id"="Inv_BID"))%>%distinct()
  
  # For some reason some pedigree names are changes in the INPUT_PEDIGREE_NAME and not found in the sql so format the Pedigree column and find it in the inb file
  
  miss1$Pedigree<- gsub("Y3_\\(","Y3_",miss1$Pedigree)
  miss1$Pedigree<- gsub("\\.)$",".",miss1$Pedigree)
  miss1$Pedigree<- gsub("%0001)$","%0001.",miss1$Pedigree)
  miss1$Pedigree<- gsub("@.0003)$","@.0003#",miss1$Pedigree)
  data4=miss1%>%left_join(inb)%>%drop_na()
  data4=data.frame(id=data4$id,GermID=data4$ProgenyGermID)
  
  if(dim(data4)[1]==0){
    alldata=alldata
  }else{
    alldata=rbind(alldata,data4)
  }
  
  miss_seas=miss_inv5%>%dplyr::filter(!id%in%alldata$id)
  miss_seas=miss_seas%>%left_join(inbredAll,by=c("id"="INPUT_BARCODE"))%>%select(id,INPUT_GERMPLASM_ID,OUTPUT_GERMPLASM_ID,INPUT_PEDIGREE_NAME)%>%distinct()
  miss3=miss_seas%>%left_join(dhSP,by=c("id"="Inv_BID"))%>%distinct()
  
  if((dim(miss3)[1] >0) == TRUE){
    write.csv(miss3,file="/mnt/DH_processing/missing_2022_parents_summer.csv",row.names = FALSE,quote = FALSE)
  }else{
    print("no missing values")
  }
  
  
  ## processing creating PCM1
  Dh_spr<-dhSP1%>%left_join(alldata,by=c(`Inv BID_Female`="id"))
  Dh_spr_1<-Dh_spr%>%left_join(alldata,by=c(`Inv BID_Male`="id"))
  
  # complete pedigree data 374 hybrids it keeps duplicated pedigree with different parent GermID for the same parent
  Dh_spr_1_1=Dh_spr_1%>%drop_na()%>%distinct()
  return(Dh_spr_1_1)
  
  
}

#########################
#### geno processing ####
#########################

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



markerdata <- s3read_using(fread, 
                           object = 'Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/MarkerMetadata.txt',
                           bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')


#load("/mnt/DH_processing/DH_2020_geno_cucumber.Rdata")
#load("/mnt/inbred_geno_cucumber.Rdata")
#load("/mnt/DH_processing/summer_DH_geno_cucumber.Rdata")#2021

### combine all data use as geno reference
#DH_2020_geno_cucumber
geno_DH <- s3read_using(fread, 
                        object = "Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/HDgenotypes.txt.gz", 
                        bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')



#summer_DH_geno_cucumber.Rdata
geno_DH_1 <- s3read_using(fread, 
                          object = "Cucumber/LongDutch/Y3/DH0_Fall2021/HDgenotypes.txt.gz", 
                          bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')


#inbred_geno_cucumber.Rdata
geno_DH_2 <- s3read_using(fread, 
                          object = "Cucumber/LongDutch/Y3/SeasonSetTrainingSet_20210825_Auto/HDgenotypes.txt.gz", 
                          bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

geno_DH=data.frame(geno_DH)
geno_DH_1=data.frame(geno_DH_1)
geno_DH_2=data.frame(geno_DH_2)

#####

geno_DH=geno_DH[!as.numeric(geno_DH$ProgenyGermID)%in%as.numeric(geno_DH_1$ProgenyGermID),]
all_geno_DH =rbind(geno_DH,geno_DH_1)

geno_DH1=all_geno_DH[!as.numeric(all_geno_DH$ProgenyGermID)%in%as.numeric(geno_DH_2$ProgenyGermID),]
all_geno_DH1 =rbind(geno_DH1,geno_DH_2)
length(unique(all_geno_DH1$ProgenyGermID))
#quality check
#length(unique(c(geno_DH$ProgenyGermID,geno_DH_1$ProgenyGermID,geno_DH_2$ProgenyGermID)))

#### function geno_proc ####
Dh_spring_1<-Dh_spr_1_1
hybridsWithBothParentsGT <- unique(Dh_spr_1_1$Pedigree)
finalFileName <- "/mnt/Y3_SummerSeasonDHSet2022_SynthHybrids.tab"

######
## This function is from the genomics team, is very slow can be optimized, parallelized
## https://stackoverflow.com/questions/14815810/append-rows-to-dataframe-using-foreach-package

geno_proc<-function(Dh_spring_1,all_geno_DH1,finalFileName){
  
  for (h in 1:nrow(Dh_spring_1)) {
    
    print(paste0("Currently working on hybrid number ",h," of ", nrow(Dh_spring_1)))
    
    
    HybridTable <- data.frame(HybridGermID=rep(hybridsWithBothParentsGT[h],nrow(markerdata)),MRN=markerdata$MRN) %>% distinct()      
    
    bothParentsImputeGT <- data.frame(MRN=markerdata$MRN)
    
    # parent1
    currentParent1 <- Dh_spring_1[h,"GermID.x"]
    currentParent1=data.frame(currentParent1)
    
    geno_DH<-data.frame(all_geno_DH1)
    fem<-geno_DH[as.numeric(geno_DH$ProgenyGermID)%in%as.numeric(currentParent1$GermID.x),]
    
    # parent2
    currentParent2 <- Dh_spring_1[h,"GermID.y"]
    tester<-geno_DH[as.numeric(geno_DH$ProgenyGermID)%in%as.numeric(currentParent2$GermID.y),]
    
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
  
}

#### function geno_format ####

geno_format<-function(file){  
  
  hybrid_geno <- read.csv(file = "/mnt/Y3_SummerSeasonDHSet2022_SynthHybrids.tab", sep="\t", header=T)
  
  hybrid_wide <- hybrid_geno %>% 
    dplyr::select(HybridGermID, MRN, Genotype) %>% 
    dcast(HybridGermID ~ MRN, value.var="Genotype")
  
  # Turn MRN to 0, 1, 2
  hybrid_wide_012 <- apply(hybrid_wide[, 2:ncol(hybrid_wide)], MARGIN = 2, 
                           FUN = function(x){
                             unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
                           })
  
  hybrid_wide_012 <- cbind(data.frame("HybridGermID"=hybrid_wide[,1]), hybrid_wide_012)
  
  # Check to see if there are any NAs in hybrid geno file
  
  tst <- hybrid_wide_012[, 3:ncol(hybrid_wide_012)]
  tst <- as.matrix(tst)
  for (i in 1:ncol(tst)){
    tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
  }
  
  tst <- as.data.frame(tst)
  hybrid_wide_012[,2:ncol(hybrid_wide_012)] <- tst
  dim(hybrid_wide_012)
  
  #s3write_using(hybrid_wide_ped, FUN = writing_csv, 
  #              object = 'Cucumber/LongDutch/Y3/SeasonSetTrainingSet_20210825_Auto/HybridGeno_wide.csv', 
  #              bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
  return(hybrid_wide_012)
}


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
  dplyr::inner_join(hybrid_wide_012, by = c('Pedigree'='HybridGermID'))

dim(geno_wide_ped)
summer_22_train<-geno_wide_ped
save(summer_22_train,file="summer_22_traincuc.Rdata")



