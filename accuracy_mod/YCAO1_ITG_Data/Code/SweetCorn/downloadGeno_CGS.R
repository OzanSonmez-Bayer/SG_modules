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

MapID <- availableTaxa[availableTaxa$description == 'Corn Concensus Map v_4',1]

map <- RetrieveMap(MapID)
head(map)
colnames(map) <- c('Chr', 'MarkerName', 'Pos')

# Translate Pedigree to Germplasm ID

source('vaultCredentials.R')
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

Crop <- 'Corn'
s3load(object = 'Corn/Corn_GermKeys.RData', bucket = "veg-apd-sdi-predictiveanalytcs-prod-reference-data")

pheno <- s3readRDS(object = 'ycao1/ITG/Corn/pheno.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
geno <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds',
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
geno_imputed <- geno
s3load(object = 'Corn/Corn_IDs.RData', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')

CropIDs_1 <- CropIDs %>% 
  select(M.GERMPLASM.X_ID, M.GERMPLASM.ORIGIN) %>% 
  unique()

origin_2020 <- read.csv('Corn/origin_2020pilot.csv')

geno_imputed_1 <- geno_imputed %>% 
  mutate(GermID = as.character(GermID)) %>% 
  left_join(CropIDs_1, by = c('GermID' = 'M.GERMPLASM.X_ID')) %>% 
  filter(as.character(M.GERMPLASM.ORIGIN) %in% as.character(origin_2020$origin))
  

# Try test retrieveCalls with GCA from PCM 0
# AWS.Authorization('ycao1')
# path <- 'ycao1/ITG_data/ProcTomato/'
# infile <- paste0('s3://','genome-analytics-perm-space/', path, 'BLUP/PGCA_PCM0_2014_2019_deep.csv')
# csvfile <- get_object(infile)
# csvfile_1 <- rawToChar(csvfile)
# con <- textConnection(csvfile_1)
# temp <- read.csv(con)
# close(con)

source('Credential/Credential_sdi.R')
# ped_list <- as.character(unique(pheno$PEDIGREE_NAME))
# Germ_IDs <- as.character(unique(geno_imputed_1$GermID))
Germ_IDs <- as.character(unique(geno_imputed$GermID))
IDs <- verifyIDs(IDList = Germ_IDs, IDType = 'GermplasmID', GoalPlatform = 'ion torrent', Crop = Crop)
batchnum <- seq(from = 1, to = length(unique(IDs[,1])), by = 100)
agg_geno <- data.frame()
LastRecord <- batchnum[1] + 99

# pull gbs data
# geno <- RetrieveCalls('ion torrent/methods/GBS/subjects/inventory/”', 
#                       CropMap = map, 
#                       GermIDList = IDs[,1])
# agg_geno <- geno[[1]]
# 
# geno_id <- merge(IDs, geno[[2]], by.x = 'GermplasmID', by.y = 'germplasm')
# colnames(agg_geno)[4:ncol(agg_geno)] <- geno_id[match(colnames(agg_geno)[4:ncol(agg_geno)], geno_id$sample), 'Pedigree']
# write.csv(agg_geno, file = 'fingerprint_data.csv', row.names = F)
# 
# IDs <- verifyIDs(IDList = ped_list, IDType = 'Pedigree', GoalPlatform = 'taqman', Crop = Crop)
# batchnum <- seq(from = 1, to = length(unique(IDs[,1])), by = 100)
# taqman_geno <- data.frame()

# LastRecord <- batchnum[1] + 99
taqman_geno <- RetrieveCalls('ion torrent/methods/GBS/subjects/inventory/',
                             CropMap = map,
                             GermIDList = unique(IDs[batchnum[1]:LastRecord,1])) 
agg_geno <- taqman_geno[[1]]
agg_id <- taqman_geno[[2]]


ptm <- proc.time()
for (i in 1238:length(batchnum)){
  print(paste0('batchnum: ', i))
  LastRecord <- batchnum[i] + 99
  if (LastRecord > length(unique(IDs[,1]))){
    LastRecord <- length(unique(IDs[,1]))
  }
  temp_taqman <- RetrieveCalls('ion torrent/methods/GBS/subjects/inventory/',
                               CropMap = map,
                               GermIDList = unique(IDs[batchnum[i]:LastRecord,1]))  
  agg_id <- rbind(agg_id, temp_taqman[[2]])
  agg_geno <- dplyr::full_join(agg_geno, temp_taqman[[1]],
                                  by = c('MarkerName' = 'MarkerName', 'Chr' = 'Chr',  'Pos' = 'Pos'))
  Sys.sleep(0.1)
}
proc.time() - ptm


# Only keep markers without NA

marker_keep <- apply(agg_geno, 
                     MARGIN = 1, 
                     FUN = function(x){
                       sum(is.na(x))
                     })
taqman_sub <- agg_geno[which(marker_keep != 10463), ]

miss_mrk <- apply(taqman_sub[,-c(1:3)], 
                  MARGIN = 2, 
                  FUN = function(x){
                    sum(is.na(x))
                  })

agg_id <- agg_id %>% 
  mutate(inventory = as.character(inventory)) %>% 
  left_join(CropIDs, by = c('inventory' = 'M.INVENTORY.X_ID'))

head(agg_id)


gbs_geno <- t(taqman_sub[, 4:ncol(taqman_sub)])

colnames(gbs_geno) <- taqman_sub$MarkerName
rownames(gbs_geno) <- colnames(taqman_sub)[-c(1:3)]

gbs_geno <- gbs_geno %>% 
  as.data.frame() %>% 
  mutate(sample = rownames(gbs_geno)) %>% 
  left_join(agg_id[, c('sample', 'M.GERMPLASM.PEDIGREE')], by = 'sample') 

gbs_geno <- gbs_geno %>% 
  select(M.GERMPLASM.PEDIGREE, sample, colnames(gbs_geno)[which(!colnames(gbs_geno) %in% c('sample', 'M.GERMPLASM.PEDIGREE'))])

gbs_geno[1:5,1:5]

s3saveRDS(gbs_geno, 
          object = 'ycao1/ITG/Corn/gbs_preimpuatation_all.rds',
          bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
