
#### function create pairwise hybrids parents
ids_par<-function(miss_inv5,inbredAll,inbredAll2,inb,dhSP,dhSP2,year,season){
  # @inbredAll sql file matched by Inv BID/Barcode
  # @inbredAll2 sql file matched by Pedigree Name
  # @inb genotyped data combined
  # @dhSP: file =  "Pedigree" "InvBID_Female" "InvBID_Male"  
  # @dhSP2: file =  "Pedigree of hybrid" "InvBID_Female" "InvBID_Male"  
  # @miss_inv5 Inv BID of parents that was input to sql
  
  Dh_fall_2<-miss_inv5%>%left_join(inbredAll,by=c("id"="INPUT_BARCODE"))
  
  ## found input germplasm
  geno1=Dh_fall_2[Dh_fall_2$INPUT_GERMPLASM_ID%in%inb$ProgenyGermID,]
  data1=data.frame(id=geno1$id,GermID=geno1$INPUT_GERMPLASM_ID)
  
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
  miss1=miss_seas%>%left_join(dhSP1,by=c("id"="Inv_BID"))%>%distinct()
  
  
  ### manual formatting step added in origin summer but not spring , should generalize
  
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
    write.csv(miss3,file= paste0("/mnt/missing_",year,"_parents_",season,".csv"),row.names = FALSE,quote = FALSE)
  }else{
    print("no missing values")
  }
  
  
  ## processing creating PCM1
  Dh_spr<-dhSP2%>%left_join(alldata,by=c(Female_InvBID="id"))
  Dh_spr_1<-Dh_spr%>%left_join(alldata,by=c(`Male tester_InvBID`="id"))
  
  # complete pedigree data 374 hybrids it keeps duplicated pedigree with different parent GermID for the same parent
  Dh_spr_1_1=Dh_spr_1%>%drop_na()%>%distinct()
  return(Dh_spr_1_1)
  
  
}


######
## This function is from the genomics team, is very slow can be optimized, parallelized
## https://stackoverflow.com/questions/14815810/append-rows-to-dataframe-using-foreach-package

geno_proc<-function(Dh_data,all_geno_DH1,finalFileName){
  #@ Dh_data list of pedigree hybrid and male and female GermID, output of ids_par
  #@ all_geno_DH1 file with all available inbreds information in an HD format
  #@ finalFileName specify under which name the file will be saved : currently uses a local folder
  
  for (h in 1:nrow(Dh_data)) {
    
    print(paste0("Currently working on hybrid number ",h," of ", nrow(Dh_data)))
    
    
    HybridTable <- data.frame(HybridGermID=rep(hybridsWithBothParentsGT[h,],nrow(markerdata)),MRN=markerdata$MRN) %>% distinct()      
    
    bothParentsImputeGT <- data.frame(MRN=markerdata$MRN)
    
    # parent1
    currentParent1 <- Dh_data[h,"GermID.x"]
    currentParent1=data.frame(currentParent1)
    
    geno_DH<-data.frame(all_geno_DH1)
    fem<-geno_DH[as.numeric(geno_DH$ProgenyGermID)%in%as.numeric(currentParent1$GermID.x),]
    
    # parent2
    currentParent2 <- Dh_data[h,"GermID.y"]
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


plan(multisession, workers = 10) 

geno_format<-function(file,crop){  
  # @file= file name to pull the HD data of hybrid same as finalFileName of geno_proc
  # @crop= crop name to pull IDs
  
  hybrid_geno <- read.csv(file = file , sep="\t", header=T)
 
  # Wide formatting of the original genotypic file
  hybrid_wide <-  hybrid_geno %>% 
    dplyr::select(HybridGermID, MRN, Genotype) %>% 
    dcast(HybridGermID ~ MRN, value.var="Genotype")
  
  FUN = function(x){
    unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
  }
  
  hybrid_wide_012<-future_apply(hybrid_wide[, 2:ncol(hybrid_wide)],MARGIN = 2, FUN=FUN)
  hybrid_wide<-data.frame(hybrid_wide)
  hybrid_wide_012<-data.frame(hybrid_wide_012)
  hybrid_wide_012 <- cbind(data.frame("ProgenyGermID"=hybrid_wide$HybridGermID), hybrid_wide_012[,1:ncol(hybrid_wide_012)])
  

  # Check to see if there are any NAs in hybrid geno file
  
  tst <- hybrid_wide_012[, 2:ncol(hybrid_wide_012)]
  tst <- as.matrix(tst)
  for (i in 1:ncol(tst)){
    tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
  }
  
  tst <- as.data.frame(tst)
  hybrid_wide_012[,2:ncol(hybrid_wide_012)] <- tst

  # Pull Reference ID from S3
  
  Crop <- crop
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
    dplyr::inner_join(hybrid_wide_012,by=c("Pedigree"="ProgenyGermID"))
  
  return(geno_wide_ped)
  
}






