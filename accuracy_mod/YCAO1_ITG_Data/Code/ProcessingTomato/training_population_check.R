# Check Genotypic data

rm(list = ls()) 
  
library(aws.s3)
library(tidyverse)
library(dplyr)
library(httr)
library(plyr)
library(data.table)


bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
crop <- 'Tomato'

# Credentials

source('vaultCredentials.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

################ Check Hybrids with Genotyped Parents #######################

years <- c('2014', '2015','2016', '2017', '2018','2019')
AllData <- data.frame()
for (selectedYear in years){
  infile <- paste0(path,crop,"/",selectedYear,"/AllYearData.RData")
  s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")
  AllData <- plyr::rbind.fill(AllData,YearPhenos)
}

tst <- read.csv('ProcTomato/PCM0_data.csv')
seasons <- c('2014:04','2015:04', '2016:04', '2017:04', '2018:04', '2019:04')

setlists <- c(unique(as.character(tst$TEST_SET_NAME)), 'PCM1', 'PCM1TI', 'PMC1LATE', 'CKRT', 'AORT', 'AO_PCM2', 'CK_PCM2')
traits_itg <- c('FBRIX', 'AVGHBRIX', 'AVJB', 'FQUAL', 'FRFRMH', 'FZUNI', 'LBRIX', 'MAT', 'MATR', 'OST', 'PLCMP', 'PLTVG', 'PYLDPA', 'SAMPUW')

rm(tst)

itg_data <- AllData %>% 
  filter(GROWSEASON %in% seasons &
           TEST_SET_NAME %in% setlists)

dim(itg_data)


itg_parent <- itg_data %>% 
  select(PEDIGREE_NAME, P1, P2) %>% 
  mutate(P2 = if_else(is.na(P2), P1, P2), 
         P1 = if_else(is.na(P1), P2, P1)) %>% 
  distinct()

dim(itg_parent)
length(unique(itg_parent$PEDIGREE_NAME))


## read hybrids ped table from S3 

hybrid_genotypePT <- s3read_using(fread, object = 'Tomato/ProcessingTomato/9Z/HybridGenotypes/Tomato_HybridPedTable_WithTrueFalseParentGT.tab',
                                    bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

head(hybrid_genotypePT)


# Read imputed gneo file fraom s3

geno <- s3read_using(fread, object = 'Tomato/ProcessingTomato/9Z/ProcessingTomato20200615-1_HDgenotypes.txt.gz', 
                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

geno_1 <- s3read_using(fread, object = 'Tomato/ProcessingTomato/9Z/imputedGeno_wide.csv',
                       bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
dim(geno_1)

# Read in refrence data to add germplasm IDs
s3load(object = 'Tomato/Tomato_IDs.RData', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')
CropIDs_sub <- CropIDs %>% 
  dplyr::select(M.GERMPLASM.PEDIGREE, M.GERMPLASM.X_ID) %>% 
  dplyr::rename(PEDIGREE_NAME = M.GERMPLASM.PEDIGREE, 
                GERMPLASM_ID = M.GERMPLASM.X_ID)

GT_germID <- as.character(unique(geno$ProgenyGermID))

GT_ped <- unique(CropIDs_sub$PEDIGREE_NAME[CropIDs_sub$GERMPLASM_ID %in% GT_germID])
length(GT_ped) # 9314

data_comparison_1 <- hybrid_genotypePT %>% 
  select(HybridPedigree, Parent1Ped, Parent2Ped, Parent1GTStatus, Parent2GTStatus, TotalParentsGT) %>% 
  full_join(itg_parent, by = c('HybridPedigree' = 'PEDIGREE_NAME')) %>% 
  mutate(P1_match = if_else(Parent1Ped == P1, 1, 0), 
         P2_match = if_else(Parent2Ped == P2, 1, 0), 
         parent_match = P1_match + P2_match,
         P1_GT = if_else(P1 %in% GT_ped, 1, 0), 
         P2_GT = if_else(P2 %in% GT_ped, 1, 0), 
         parent_GT = P1_GT + P2_GT,
         P1_check = if_else(Parent1Ped %in% GT_ped, 1, 0), 
         P2_check = if_else(Parent2Ped %in% GT_ped, 1, 0), 
         parent_check = P1_check + P2_check)


dim(data_comparison_1)
head(data_comparison_1)

table(data_comparison_1$parent_GT)

# 0    1    2 
# 132 1291 6399

table(data_comparison_1$TotalParentsGT)

# 0    1    2 
# 68 1273 6392

table(data_comparison_1$parent_match) # Parents provided by genomics team are not correct 

table(data_comparison_1$parent_check)


parents_list <- unique(c(itg_parent$P1, itg_parent$P2))
parents_list_genomics <- unique(c(data_comparison_1$Parent1Ped, data_comparison_1$Parent2Ped))

length(intersect(parents_list, GT_ped))
length(intersect(parents_list_genomics, GT_ped))


######################### Check Genotypes by Traits #####################

trait_df <- itg_data %>% 
  filter(OBSRVTN_REF_CD %in% traits_itg) %>% 
  select(PEDIGREE_NAME, P1, P2, OBSRVTN_REF_CD) %>% 
  mutate(P2 = if_else(is.na(P2), P1, P2)) %>% 
  distinct()

dim(trait_df) # 510148 4


trait_df_1 <- trait_df %>%
  dplyr::mutate(P1_genotyped = if_else(as.character(P1) %in% GT_ped, 1, 0),
                P2_genotyped = if_else(as.character(P2) %in% GT_ped, 1, 0),
                parents_genotyped = P1_genotyped + P2_genotyped)

table(trait_df_1$OBSRVTN_REF_CD, trait_df_1$parents_genotyped)

ct_trait <- as.data.frame(table(trait_df_1$OBSRVTN_REF_CD, trait_df_1$parents_genotyped)) %>% 
  reshape2::dcast(Var1 ~ Var2,value.var = 'Freq' ) %>% 
  dplyr::rename(Trait = Var1)

ct_trait
