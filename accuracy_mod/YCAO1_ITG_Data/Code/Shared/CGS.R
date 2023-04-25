
################################ CGS FUNCTIONS ########################

library(httr)
library(jsonlite)
library(ggplot2)
library(plotly)
library(doBy)
library(hashmap)

# AvailablePlatforms <- hashmap( keys = c('aggregated','ion torrent','infinium'),
#                                values = c('aggregated/methods/consolidation','ion torrent/methods/GBS','infinium/methods/component')
# ) # END of Platform table


### FROM AUTHENTICATOR 0.2
#########################################
createPingToken <- function (passwordFlow = FALSE, maxRetries = 3) 
{
  credentials = Sys.getenv("PATH_TO_CREDENTIALS", unset = NA)
  if (Sys.getenv("PATH_TO_CREDENTIALS") != "") {
    credentials <- read.table(Sys.getenv("PATH_TO_CREDENTIALS"), 
                              header = FALSE)
    client_secret <- as.character(credentials[1, ])
    client_id <- as.character(credentials[2, ])
    url <- as.character(credentials[3, ])
    if (passwordFlow) {
      username <- as.character(credentials[4, ])
      password <- as.character(credentials[5, ])
      if (is.na(username) || username == "" || is.na(password) || 
          password == "") {
        warning("password-grant flow requires username and password")
      }
      if (is.na(username) || username == "") {
        warning("username is not provided in your credentials file")
        username = ""
      }
      if (is.na(password) || password == "") {
        warning("password is not provided in your credentials file")
        password = ""
      }
    }
  }
  else {
    client_id <- Sys.getenv("CLIENT_ID")
    client_secret <- Sys.getenv("CLIENT_SECRET")
    url <- Sys.getenv("PING_URL")
    if (passwordFlow) {
      username <- Sys.getenv("USERNAME")
      password <- Sys.getenv("PASSWORD")
      if (is.na(username) || username == "" || is.na(password) || 
          password == "") {
        warning("password-grant flow requires username and password")
      }
      if (is.na(username) || username == "") {
        warning("Environmental variable USERNAME is not set")
        username = ""
      }
      if (is.na(password) || password == "") {
        warning("Environmental variable PASSWORD is not set")
        password = ""
      }
    }
  }
  auth_header <- paste0("Basic ", base64enc::base64encode(charToRaw(paste0(client_id, 
                                                                           ":", client_secret))))
  attempts <- 0
  repeat {
    if (passwordFlow) {
      response <- httr::POST(url, httr::add_headers(Authorization = auth_header), 
                             httr::add_headers(`Accept-Encoding` = "gzip, deflate"), 
                             httr::user_agent("authenticator"), httr::content_type("application/x-www-form-urlencoded;charset=UTF-8"), 
                             body = list(grant_type = "password", username = username, 
                                         password = password), encode = "form")
    }
    else {
      response <- httr::POST(url, httr::add_headers(Authorization = auth_header), 
                             httr::add_headers(`Accept-Encoding` = "gzip, deflate"), 
                             httr::user_agent("authenticator"), httr::content_type("application/x-www-form-urlencoded;charset=UTF-8"), 
                             body = list(grant_type = "client_credentials"), 
                             encode = "form")
    }
    attempts <- attempts + 1
    if (httr::status_code(response) == 200) {
      break
    }
    if (attempts >= maxRetries) {
      stop("Maximum number of retries to receive OAuth token exceed")
    }
  }
  BEARER_TOKEN <<- paste0("Bearer ", httr::content(response)$access_token)
  return(BEARER_TOKEN)
}

#createPingToken

pingSecuredGet <- function(url, userAgent = "authenticator", query = NULL, timeOut = 60) 
{
  url <- URLencode(url)
  response <- httr::GET(url, httr::add_headers(Authorization = BEARER_TOKEN), 
                        httr::user_agent(userAgent), httr::timeout(timeOut), 
                        query = query)
  if (response$status_code %in% c(302, 401, 403)) {
    httr::handle_reset(url)
    createPingToken()
    response <- pingSecuredGet(url, userAgent, timeOut)
  }
  return(response)
}

######################################### CGS REQUESTS #########################################

###### API ACCESS

#### RETRIEVING MAP
RetrieveMap <- function (MapID) {
  
  #MapID <- AvailableTaxaTable[grep(InputCrop, AvailableTaxaTable$name, ignore.case = TRUE),1]
  
  MapCallString <- paste0("https://genetics.ag/genetic-maps/v1/",MapID)
  
  createPingToken()
  MapJson <- pingSecuredGet(MapCallString)
  MapContigNames <- as.numeric(fromJSON(content(MapJson, as = "text", encoding = "UTF-8"))$contigMaps$contigName)
  MapContigPositions <- fromJSON(content(MapJson, as = "text", encoding = "UTF-8"))$contigMaps$positions
  
  ## Concatenating all the list elements into a single object
  CompleteMap <- data.frame()
  for (ContigIndex in 1:length(MapContigNames)) {
    CompleteMap <- rbind(CompleteMap,cbind(rep(MapContigNames[ContigIndex], nrow(MapContigPositions[[ContigIndex]])), MapContigPositions[[ContigIndex]]))
  } # End getting all markers
  
  colnames(CompleteMap) <- c('Chr','Pos', 'MarkerName')
  
  #CompleteMap$Chr <- as.numeric(CompleteMap$Chr)
  
  CompleteMap[with(CompleteMap, order(Chr,Pos)),]
  #CompleteMap
  
}


#### GET MAP STATS
GetMapStats <- function(Map){
  
  CountPerChr <- table(Map$Chr)
  
  chr_len = sapply(split(Map$Pos, Map$Chr),max)
  
  Keys <- c('Number of chromosomes', 'Total number of markers', 'Total size in cM', 'Min number of markers per chromosome', 'Max number of markers per chromosome')
  Values <- c(length(CountPerChr), nrow(Map), sum(chr_len), min(CountPerChr), max(CountPerChr))
  
  data.frame(Keys,Values)
  
}


#### RETRIEVING CALLS
RetrieveCalls <- function(Platform,CropMap,GermIDList) {
  #RetrieveCalls <- function(Map,GermIDList) {
  
  CallTable <- data.frame(CropMap)
  #CallTable <- data.frame()
  
  SubjectTable <- data.frame()
  
  
    
  for(GermID in GermIDList) {
      
    #print(GermID)
      
      
    CallString = paste0("https://genetics.ag/marker-call-sets/v1/platforms/",Platform, GermID,"/call-sets")
    #CallString = paste0("https://genetics.ag/marker-call-sets/v1/platforms/aggregated/methods/consolidation/subjects/germplasm/", GermID,"/call-sets")
      
    #CallTable <- rbind(CallTable,c(GermID,CallString))
    MarkerSetJson <- pingSecuredGet(CallString, timeOut = 40)
      
    if(MarkerSetJson$status_code == '200') {
      MarkerSet <- fromJSON(content(MarkerSetJson, as = "text", encoding = "UTF-8"))
        
      RetrievedSubject <- data.frame(matrix(unlist(MarkerSet$subjects), nrow=length(MarkerSet$subjects),byrow=T))[,c(3,4)]
      colnames(RetrievedSubject) <- MarkerSet$subjects[[1]]$subjectType
        
      SubjectTable <- rbind(SubjectTable, RetrievedSubject)
        
      for(i in 1:length(MarkerSet$subjects)) {
          
        SampleID <- MarkerSet$subjects[[i]][which(MarkerSet$subjects[[i]]$subjectType == 'sample'),'subjectId']
          
        MarkerCalls <- data.frame(
          MarkerName = as.data.frame(MarkerSet$calls[[i]])[ , "markerName"],
          Calls = unlist(lapply(as.data.frame(MarkerSet$calls[[i]])[ , "alleleCalls"], function(x) paste(unlist(x), collapse = "|"))),
          stringsAsFactors = FALSE)
          
          
          
        colnames(MarkerCalls)[2] <- SampleID
          
        CallTable <- merge(x=CallTable, y=MarkerCalls, by=0,by.x="MarkerName",by.y="MarkerName", all = TRUE)
          
          #CallTable <- cbind(CallTable, MarkerCalls$Calls[match(CallTable$MarkerName,MarkerCalls$MarkerName)])
          
          # if(length(which(is.na(match(MarkerCalls$MarkerName, CallTable$MarkerName))))>0){
          #
          #   MissingMarkers <- MarkerCalls[which(is.na(match(MarkerCalls$MarkerName, CallTable$MarkerName))),]
          #
          #   NumMissing <- nrow(MissingMarkers)
          #
          #   if(ncol(CallTable)==4) {
          #     MissingDataset <- cbind(rep(0,NumMissing),rep(0,NumMissing), MissingMarkers$MarkerName, MissingMarkers$Calls)
          #   } else {
          #     MissingDataset <- cbind(rep(0,NumMissing),rep(0,NumMissing), MissingMarkers$MarkerName, matrix(rep(NA, NumMissing*(ncol(CallTable)-4)), ncol = ncol(CallTable)-4), MissingMarkers[,SampleID])
          #   }
          #
          #   colnames(MissingDataset) <- colnames(CallTable)
          #
          #   CallTable <- rbind(CallTable, MissingDataset)
          #
          # } #end if there are missing markers
          
      } #end for looping across subjects
        
    } else {
        
        #CallTable <- cbind(CallTable, rep(NA,nrow(CallTable)))
        
    } # End if call worked
      
  } # End of loop through accessions
  
  #colnames(CallTable) <- c("Chr","Pos","MarkerName", GermIDList)
  
 list(CallTable, SubjectTable)
  
  #CountMissingPerMarkers <- apply(CallTable, 1, function(x) length(which(is.na(x))))
  #CallTable[which(CountMissingPerMarkers<length(GermIDList)),]
  
}


#### CALCULATE MAF

CalculateMAF <- function(GenoCalls) {
  
  #GenoCalls <- c("CC","GG",NA,"CC","CC","GG","CC",NA,"GG","GG","GG",NA,"GG")
  #GenoCalls <- c("CC","GG","CC")
  
  if(length(which(is.na(GenoCalls)))>0){
    GenoCalls <- GenoCalls[-(which(is.na(GenoCalls)))]
  } #if data has NA's
  
  Alleles <- sort(table(GenoCalls), decreasing = FALSE)
  
  if(length(Alleles)==1) {
    
    0
    
  } else {
    
    Alleles[1]/length(GenoCalls)
    
  } # End if there's only one allele
  
} # END MAF Calculation



######################################### ID TRANSLATION #########################################


verifyIDs <- function(IDList, IDType, GoalPlatform, Crop) {
  
  translatedIDList <- data.frame()
  
  UniqueIDs <- unique(IDList)
  
  IDColnames <- list(
    "InventoryID" = "M.INVENTORY.X_ID",
    "InventoryBID" = "M.INVENTORY.X_BID",
    "Pedigree" = "M.GERMPLASM.PEDIGREE",
    "GermplasmID" = "M.GERMPLASM.X_ID"
  )
  
    # source('Credential/Credential_sdi.R')
    
    ############ FOR AGGREGATED DATA ############
    if(GoalPlatform == 'aggregated')
    {
      
      if(IDType == 'GermplasmID') {
        
        translatedIDList <- data.frame(M.GERMPLASM.X_ID = UniqueIDs)
        colnames(translatedIDList) <- c("GermplasmID")
        
      } else if (IDType %in% c('Pedigree')) {
        
        GermKeysObject <- paste(Crop, "/",Crop, "_GermKeys.RData", sep="")
        s3load(object = GermKeysObject, bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")
        
        print(GermKeysObject)
        
        translatedIDList <- CropGermKeys[which(CropGermKeys$M.GERMPLASM.PEDIGREE %in% UniqueIDs),c("M.GERMPLASM.X_ID","M.GERMPLASM.PEDIGREE")]
        colnames(translatedIDList) <- c("GermplasmID","Pedigree")
        
      } else if (IDType %in% c('InventoryID', 'InventoryBID')) {
        
        IDKeysObject <- paste(Crop, "/",Crop, "_IDs.RData", sep="")
        s3load(object = IDKeysObject, bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")
        
        ColToSearch <- which(colnames(CropIDs) == IDColnames[[IDType]])
        
        translatedIDList <- CropGermKeys[which(CropGermKeys[,ColToSearch] %in% UniqueIDs),c("M.GERMPLASM.X_ID",IDColnames[[IDType]])]
        colnames(translatedIDList) <- c("GermplasmID",IDType)
        
      } #End of checking provided ID
      
      
      ############ FOR COMPONENT DATA ############
    } else if (GoalPlatform %in% c('ion torrent','infinium', 'taqman'))
    {
      
      if(IDType == 'InventoryID') {
        
        translatedIDList <- data.frame(M.INVENTORY.X_ID = UniqueIDs)
        colnames(translatedIDList) <- c("InventoryID")
        
      } else if (IDType %in% c('GermplasmID','Pedigree', 'InventoryBID')) {
        
        IDKeysObject <- paste(Crop, "/",Crop, "_IDs.RData", sep="")
        s3load(object = IDKeysObject, bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")
        
        ColToSearch <- which(colnames(CropIDs) == IDColnames[[IDType]])
        
        translatedIDList <- CropIDs[which(CropIDs[,ColToSearch] %in% UniqueIDs),c("M.INVENTORY.X_ID",IDColnames[[IDType]])]
        colnames(translatedIDList) <- c("InventoryID",IDType)
        
        
        if(length(which(translatedIDList$M.INVENTORY.X_ID == "")) > 0 ){
          translatedIDList <- translatedIDList[-(which(translatedIDList$M.INVENTORY.X_ID == "")),]
        }
        
        #translatedIDList <- translatedIDList[-(which(is.na(translatedIDList$M.INVENTORY.X_ID))),]
        
      } #End of checking provided ID
      
    } else {
      
      print("Error: Goal platform not recognized")
      
    } ## End defining platform type

  
  translatedIDList
  
} #End of Verify IDs function
