library(aws.s3)

source('Code/Shared/CGS.R')
source('Code/Shared/misFunctions.R')
source('Credential/Credential_sdi.R')
Sys.setenv(CLIENT_SECRET = FINGERPRINT_CLIENT_SECRET)
Sys.setenv(CLIENT_ID = FINGERPRINT_CLIENT_ID)
Sys.setenv(PING_URL = FINGERPRINT_PING_URL)

createPingToken()

availableTaxaJson <- pingSecuredGet(URLencode('https://genetics.ag/genetic-maps/v1/'))

availableTaxa <- fromJSON(content(availableTaxaJson, as = "text", encoding = "UTF-8"))

MapID <- availableTaxa[availableTaxa$description == 'Tomato Consensus Map, June 2014',1]

map <- RetrieveMap(MapID)
head(map)
colnames(map) <- c('Chr', 'MarkerName', 'Pos')

# Translate Pedigree to Germplasm ID
Crop <- 'Tomato'
GermKeysObject <- paste("/shared/Veg_Phenotypes/Ref_Data/",Crop, "/",Crop, "_GermKeys.RData", sep="")
opts <- list()
opts$headers = c('x-amz-server-side-encryption' = 'AES256')
s3load(object = GermKeysObject, bucket = "genome-analytics-perm-space", opts = opts)


# Try test retrieveCalls with GCA from PCM 0
AWS.Authorization('ycao1')
path <- 'ycao1/ITG_data/ProcTomato/'
infile <- paste0('s3://','genome-analytics-perm-space/', path, 'BLUP/PGCA_PCM0_2014_2019_deep.csv')
csvfile <- get_object(infile)
csvfile_1 <- rawToChar(csvfile)
con <- textConnection(csvfile_1)
temp <- read.csv(con)
close(con)

source('Credential/Credential_sdi.R')
ped_list <- as.character(unique(temp$PEDIGREE_NAME))
IDs <- verifyIDs(IDList = ped_list, IDType = 'Pedigree', GoalPlatform = 'aggregated', Crop = Crop)
batchnum <- seq(from = 1, to = length(unique(IDs[,1])), by = 100)
agg_geno <- data.frame()
LastRecord <- batchnum[1] + 99


geno <- RetrieveCalls('aggregated/methods/consolidation/subjects/germplasm/', 
                     CropMap = map, 
                     GermIDList = unique(IDs[,1]))
agg_geno <- geno[[1]]

geno_id <- merge(IDs, geno[[2]], by.x = 'GermplasmID', by.y = 'germplasm')
colnames(agg_geno)[4:ncol(agg_geno)] <- geno_id[match(colnames(agg_geno)[4:ncol(agg_geno)], geno_id$sample), 'Pedigree']
write.csv(agg_geno, file = 'fingerprint_data.csv', row.names = F)

IDs <- verifyIDs(IDList = ped_list, IDType = 'Pedigree', GoalPlatform = 'taqman', Crop = Crop)
batchnum <- seq(from = 1, to = length(unique(IDs[,1])), by = 100)
taqman_geno <- data.frame()

LastRecord <- batchnum[1] + 99
taqman_geno <- RetrieveCalls('taqman/methods/taqman/subjects/inventory/',
                             CropMap = map,
                             GermIDList = unique(IDs[batchnum[1]:LastRecord,1])) 
taqman


ptm <- proc.time()
for (i in 1:length(batchnum)){
  LastRecord <- batchnum[i] + 99
  if (LastRecord > length(unique(IDs[,1]))){
    LastRecord <- length(unique(IDs[,1]))
  }
    temp_taqman <- RetrieveCalls('taqman/methods/taqman/subjects/inventory/',
                                 CropMap = map,
                                 GermIDList = unique(IDs[batchnum[i]:LastRecord,1]))  
    taqman_geno <- rbind(taqman_geno, temp_taqman[[2]])
    Sys.sleep(0.1)
}
proc.time() - ptm


# Import a list of inventory IDs 
inventory_df <- read_excel('C:/Users/YCAO1/Desktop/P-1/ProcTomato/9Z_TaqMan_Mabmart.xlsx', sheet = '9Z_TaqMan_Mabmart')

inventory_id <- inventory_df$INDIVIDUAL_ID...1
length(inventory_id)
str(inventory_id)

# Pull taqman data based on sample id
taqman_geno <- RetrieveCalls(Platform = 'taqman/methods/taqman/subjects/sample/', 
                              CropMap = map, 
                              GermIDList = inventory_id[1])


batchnum <- seq(from = 1, to = length(inventory_id), by = 100)
taqman_geno <- data.frame()
LastRecord <- batchnum[1] + 99
taqman_geno <- RetrieveCalls('taqman/methods/taqman/subjects/sample/',
                             CropMap = map,
                             GermIDList = unique(inventory_id[batchnum[1]:LastRecord])) 
taqman_geno <- taqman_geno[[1]]

ptm <- proc.time()
for (i in 2:length(batchnum)){
  print(paste("batch number: ", i))
  LastRecord <- batchnum[i] + 99
  if (LastRecord > length(unique(inventory_id))){
    LastRecord <- length(unique(inventory_id))
  }
    temp_taqman <- RetrieveCalls('taqman/methods/taqman/subjects/sample/',
                                 CropMap = map,
                                 GermIDList = unique(inventory_id[batchnum[i]:LastRecord]))  
    taqman_geno <- cbind(taqman_geno, temp_taqman[[1]][, -c(1:3)])
    Sys.sleep(0.1)
}
proc.time() - ptm


# Only keep markers without NA

marker_keep <- apply(taqman_geno, 
                     MARGIN = 1, 
                     FUN = function(x){
                       sum(is.na(x))
                     })
taqman_sub <- taqman_geno[which(marker_keep != 7560), ]

miss_mrk <- apply(taqman_sub[,-c(1:3)], 
                     MARGIN = 2, 
                     FUN = function(x){
                       sum(is.na(x))
                     })


#################### Inventory ID 
inventory_dat <- read_excel('C:/Users/YCAO1/Desktop/P-1/ProcTomato/9Z_Mabmart_TaqManIDPulls.xlsx')
inventory_dat <- inventory_dat %>% 
  filter(`InvBIDs Sent to TQ` == TRUE)

inventory_ls <- inventory_dat$INVENTORY_ID

batchnum <- seq(from = 1, to = length(inventory_ls), by = 100)
taqman_geno <- data.frame()

LastRecord <- batchnum[1]+ 99
taqman_geno <- RetrieveCalls('taqman/methods/taqman/subjects/inventory/',
                             CropMap = map,
                             GermIDList = unique(inventory_ls[batchnum[1]:LastRecord])) 
taqman_geno <- taqman_geno[[1]]
# save(taqman_geno, file = 'C:/Users/YCAO1/Desktop/P-1/ProcTomato/taqman_geno_1.RData')

ptm <- proc.time()
for (i in 2:length(batchnum)){
  print(paste('batch number: ', i))
  LastRecord <- batchnum[i] + 99
  if (LastRecord > length(unique(inventory_ls))){
    LastRecord <- length(unique(inventory_ls))
  }
  temp_taqman <- RetrieveCalls('taqman/methods/taqman/subjects/inventory/',
                               CropMap = map,
                               GermIDList = unique(inventory_ls[batchnum[i]:LastRecord]))  
  taqman_geno <- dplyr::full_join(taqman_geno, temp_taqman[[1]],
                                  by = c('MarkerName' = 'MarkerName', 'Chr' = 'Chr',  'Pos' = 'Pos'))
  # taqman_geno <- temp_taqman[[1]]
  # save(taqman_geno, 
  #      file = paste0(paste('C:/Users/YCAO1/Desktop/P-1/ProcTomato/taqman_geno', 1, sep = "_"), '.RData', collapse = ''))
  # rm(taqman_geno)
  rm(temp_taqman)
  Sys.sleep(0.1)
}
proc.time() - ptm



# Combine all the taqman data together

setwd('C:/Users/YCAO1/Desktop/P-1/ProcTomato')
file_ls <- list.files(path = 'C:/Users/YCAO1/Desktop/P-1/ProcTomato', pattern= '.RData')

load("taqman_geno.RData")
# duplicated columns
sample_ID <- unlist(lapply(colnames(taqman_geno), function(x){unlist(strsplit(x, '\\.'))[1]}))
# sample_ID <- unique(sample_ID)
colnames(taqman_geno) <- sample_ID
taqman_geno <- taqman_geno[, !duplicated(colnames(taqman_geno))]
procTomato_taqman <- taqman_geno

for (file_name in 2:length(file_ls)){
  print(file_name)
  load(file_ls[file_name])
  procTomato_taqman <- dplyr::full_join(procTomato_taqman, taqman_geno, 
                                        by = c('MarkerName' = 'MarkerName', 'Chr' = 'Chr',  'Pos' = 'Pos'))
  rm(taqman_geno)
}


sample_ID <- unlist(lapply(colnames(procTomato_taqman), function(x){unlist(strsplit(x, '\\.'))[1]}))
colnames(procTomato_taqman) <- sample_ID
# sample_ID <- unique(sample_ID)

procTomato_taqman <- procTomato_taqman[, !duplicated(colnames(procTomato_taqman))]


marker_keep <- apply(procTomato_taqman, 
                     MARGIN = 1, 
                     FUN = function(x){
                       sum(is.na(x))
                     })
taqman_sub <- procTomato_taqman[which(marker_keep != 21173), ]

miss_mrk <- apply(taqman_sub, 
                  MARGIN = 2, 
                  FUN = function(x){
                    sum(is.na(x))
                  })

save(procTomato_taqman, file = 'taqman.RData')


s3write_using(taqman_sub, FUN = save, 
              bucket = 'genome-analytics-perm-space', 
              object = 'ycao1/ITG_data/ProcTomato/Data/Geno/Taqman.RData')
