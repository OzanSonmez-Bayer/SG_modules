rm(list = ls())

library(aws.s3)
library(tidyverse)
library(dplyr)
library(httr)
library(plyr)
library(data.table)


bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
crop <- 'Corn'

# Credentials

source('vaultCredentials.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")



######################## Check the number of hybrids with genotyped parents ##################

# Getting Phenotypic Data
X6S_TestSets <- readxl::read_excel('/rstudio-efs/home/ycao1/ITG_Data/Corn/6S_TestSets.xlsx')
setlists <- unique(X6S_TestSets$TestSetName)

# P1 and P2 are missing from 2014
years <- c('2014', '2015','2016', '2017', '2018','2019')
AllData <- data.frame()
for (selectedYear in years){
  infile <- paste0(path,crop,"/",selectedYear,"/AllYearData.RData")
  s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")
  AllData <- plyr::rbind.fill(AllData,YearPhenos)
}

seasons <- c('2014:04','2015:04', '2016:04', '2017:04', '2018:04', '2019:04')

# traits_itg <- c('HSC', 'RTLR', 'SC_TF', 'SDV', 'ELEN', 'EDIA', 'SC_SE', 'SC_HE', 'QUAL', 'S50D', 'SC_PR')
traits_itg <- c('SC_TF', 'EDIA', 'ELEN', 'HSC', 'PA', 'SC_HE', 'SC_SE', 'SC_PR', 'SDV', 'ESHK', 'TILLR',
                'RTLR', 'QUAL', 'SC_FG', 'SC_FR', 'S50', 'S50_BE', 'LDSR_SETOTU', 'LDSR_PUCCSO', 'PCTRC')

# traits_itg <- 'PCTRC'

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



hybrid_genotypePT <- readxl::read_excel('Corn/hybrid_genotypePT.xlsx')
genotyped_parents <- unique(c(hybrid_genotypePT$Parent1Ped[hybrid_genotypePT$Parent1GTStatus == 1], 
                              hybrid_genotypePT$Parent2Ped[hybrid_genotypePT$Parent2GTStatus == 1]))  
length(genotyped_parents)

genotyped_parents_ID <- unique(c(hybrid_genotypePT$Parent1[hybrid_genotypePT$Parent1GTStatus == 1], 
                                 hybrid_genotypePT$Parent2[hybrid_genotypePT$Parent2GTStatus == 1]))  

length(genotyped_parents_ID)

hybridGT <- s3read_using(readxl::read_xlsx, bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data',
                         object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/6S_HybridInvPedTableWithPeds.xlsx')

geno_1 <- s3read_using(fread, bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data', 
                       object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Sweetcorn20200703-1_HDgenotypes.txt.gz')



# Try to match parents with germplasm ID

s3load(object = 'Corn/Corn_IDs.RData', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')
CropIDs_sub <- CropIDs %>% 
  dplyr::select(M.GERMPLASM.PEDIGREE, M.GERMPLASM.X_ID) %>% 
  dplyr::rename(PEDIGREE_NAME = M.GERMPLASM.PEDIGREE, 
                GERMPLASM_ID = M.GERMPLASM.X_ID)


GT_germID <- as.character(unique(geno_1$ProgenyGermID))

GT_ped <- unique(CropIDs_sub$PEDIGREE_NAME[CropIDs_sub$GERMPLASM_ID %in% GT_germID])
length(GT_ped) # 9314

data_comparison_1 <- hybrid_genotypePT %>% 
  select(HybridPed, Parent1Ped, Parent2Ped, Parent1GTStatus, Parent2GTStatus, TotalParentsGT) %>% 
  full_join(itg_parent, by = c('HybridPed' = 'PEDIGREE_NAME')) %>% 
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

table(data_comparison_1$parent_GT) 
# 0     1     2 
# 138 10590 14694 

table(data_comparison_1$TotalParentsGT)
# 0     1     2 
# 89 10677 14539

table(data_comparison_1$parent_match)
# 0     1     2 
# 1     5 25296 

table(data_comparison_1$parent_check)
# 0     1     2 
# 206 10576 14640 


# Non-overlapped hybrids

hybrids_diff <- setdiff(itg_parent$PEDIGREE_NAME, unique(hybrid_genotypePT$HybridPed))

length(itg_parent$PEDIGREE_NAME) # 24986
length(unique(hybrid_genotypePT$HybridPed)) # 24868

length(hybrids_diff)

sum(unique(hybrid_genotypePT$HybridPed) %in% itg_parent$PEDIGREE_NAME) # itg_parent has all the hybrids 

# non-overlap parents
parents_list_itg <- unique(c(itg_parent$P1, itg_parent$P2))
parents_list_itg <- parents_list_itg[!is.na(parents_list_itg)]

parents_list_genomic <- unique(c(hybrid_genotypePT$Parent1Ped, hybrid_genotypePT$Parent2Ped))
parents_list_genomic <- parents_list_genomic[!is.na(parents_list_genomic)]

length(parents_list_itg) # 15117
length(genotyped_parents) # 9157
length(parents_list_genomic) # 15099

parent_diff <- setdiff(parents_list_itg, parents_list_genomic)
length(parent_diff)

sum(parents_list_genomic %in% parents_list_itg) # All the parents provided by Genomics are in the parents_list_itg 

# Parents expected to be genotyped

expected_parent_GT <- setdiff(parents_list_itg, genotyped_parents)
length(expected_parent_GT) # 5960


xlsx::write.xlsx(as.data.frame(table(data_comparison_1$parent_GT) ), 
                 file = 'Corn/training_population_comparison.xlsx', 
                 sheetName = 'overview', append = F,row.names = F)
xlsx::write.xlsx(hybrids_diff, file = 'Corn/training_population_comparison.xlsx', sheetName = 'non-overlap hybrids', append = T,row.names = F)
xlsx::write.xlsx(parent_diff, file = 'Corn/training_population_comparison.xlsx', sheetName = 'non-overlap parents', append = T,row.names = F)
xlsx::write.xlsx(expected_parent_GT, file = 'Corn/training_population_comparison.xlsx', sheetName = 'Parent not genotyped', append = T,row.names = F)


################ Check GT parents by trait ##############

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
  rename(Trait = Var1)

ct_trait

write.csv(ct_trait, file = 'Corn/GenotypedParents_by_trait.csv', row.names = F)


################ Add Germplasm IDs ###########

germ_df <- itg_parent %>% 
  left_join(CropIDs_sub, by = 'PEDIGREE_NAME') %>% 
  rename(HybridGermID = GERMPLASM_ID) %>% 
  distinct() %>% 
  left_join(CropIDs_sub, by = c('P1' = 'PEDIGREE_NAME')) %>% 
  rename(Parent1GermID = GERMPLASM_ID) %>% 
  distinct() %>% 
  left_join(CropIDs_sub, by = c('P2' = 'PEDIGREE_NAME')) %>% 
  rename(Parent2GermID = GERMPLASM_ID) %>% 
  distinct()

germ_df <- germ_df %>% 
  mutate(HybridGermID = as.character(HybridGermID),
         Parent1GermID = as.character(Parent1GermID),
         Parent2GermID = as.character(Parent2GermID))

dim(germ_df)
length(unique(germ_df$PEDIGREE_NAME))
length(unique(germ_df$HybridGermID))

unique_hybrid <- unique(germ_df$PEDIGREE_NAME)
unique_parents <- unique(c(germ_df$P1, germ_df$P2))

xlsx::write.xlsx(germ_df, file = 'Corn/Hybrid_parent.xlsx', append = F, row.names = F, sheetName = 'hybrid_parent')
xlsx::write.xlsx(unique_hybrid, file = 'Corn/Hybrid_parent.xlsx', append = T, row.names = F, sheetName = 'hybrid')
xlsx::write.xlsx(unique_parents, file = 'Corn/Hybrid_parent.xlsx', append = T, row.names = F, sheetName = 'parent')
