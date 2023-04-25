# Augemented RRBLUP: snps with large effects are fixed effects. 

library(aws.s3)
library(rrBLUP)
library(dplyr)
library(httr)
library(tidyverse)
library(data.table)


source('vaultCredentials.R')

bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
crop <- 'Corn'

# Credentials

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")


########################################### Import pheno and geno ###########################

pheno <- s3readRDS(object = 'ycao1/ITG/Corn/pheno.rds',
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

# pheno from training population 
pheno <- pheno %>% 
  filter(GROWSEASON != '2020:04')

geno <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds', 
                  bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
# Import marker metadata

marker_meta <- s3read_using(fread, object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Archive/Sweetcorn20200703-1_MarkerMetadata.txt', 
                            bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

# Test wtih one trait

dat <- pheno %>% 
  filter(OBSRVTN_REF_CD == 'SDV') %>% 
  select(PEDIGREE_NAME, TRAIT_VALUE) %>% 
  filter(PEDIGREE_NAME != '') %>%
  mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
  group_by(PEDIGREE_NAME) %>% 
  summarise(lsmean = mean(TRAIT_VALUE, na.rm = T))

geno_sub <- geno %>% 
  filter(Pedigree %in% as.character(dat$PEDIGREE_NAME))

geno_sub <- geno_sub[!duplicated(geno_sub$Pedigree),]

geno_sub[, 3:ncol(geno_sub)] <- geno_sub[, 3:ncol(geno_sub)] - 1

geno_sub[1:5,1:5]


################## Get Trait linked Markers ###################

# s3load(object = 'GPC_Data/Corn/Corn_GPC.RData',
#        bucket = 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data')

################################################################

###### Calculate Missing Markers 

geno_imputed <- A.mat(geno_sub[, 3:ncol(geno_sub)], return.imputed = T)
geno_no_missing <- geno_imputed$imputed

geno_df <- data.frame(Pedigree = geno_sub$Pedigree, geno_no_missing)

geno_pheno <- dat %>% 
  inner_join(geno_df, by = c('PEDIGREE_NAME' = 'Pedigree'))

geno_no_missing <- geno_pheno[, 3:ncol(geno_pheno)]

# reformat geno for rrBLUP

geno_t <- as.data.frame(t(geno_no_missing))
colnames(geno_t) <- geno_pheno$PEDIGREE_NAME
geno_t$MRN <- rownames(geno_t)


# SNP positions

geno_pos <- marker_meta %>%
  select(MRN, LG, cM_times_1e6)



# geno_cln <- data.frame(MRN = as.character(rownames(geno_sub_t)))
geno_cln <- geno_pos %>%
  inner_join(geno_t, by = 'MRN') %>%
  rename(pos = cM_times_1e6,
         Chr = LG)


# geno_pheno <- data.frame(geno_pheno[, c('PEDIGREE_NAME', 'lsmean')], geno_cln)


# 
# geno_cln <- cbind(geno_cln, geno_sub_t)
# 
# # Match pedigrees in pheno with peds in geno
# 
# pheno_sub <- dat %>% 
#   filter(PEDIGREE_NAME %in% colnames(geno_cln)[4:ncol(geno_cln)])


# rm data

rm(geno_df)
rm(geno_no_missing)
rm(geno_pos)
rm(geno_sub)
rm(geno_t)
gc()



######################################### GWAS with RRBLUP ############################

set.seed(2020)


# Create 5-fold Cross Validation 

library(caret)

flds <- createFolds(geno_pheno$PEDIGREE_NAME, k = 5, list = TRUE, returnTrain = TRUE)

for (i_fld in 1:length(flds)){
  pheno_df <- geno_pheno[flds[[i_fld]], ] %>% 
    select(PEDIGREE_NAME, lsmean)
  
  pheno_df <- data.frame(pheno_df)
  
  colnames_selected <- colnames(geno_cln)[which(colnames(geno_cln) %in% c('MRN', 'Chr', 'pos', pheno_df$PEDIGREE_NAME))]
  geno_df <- geno_cln %>% 
    select(colnames_selected)
  
  system.time(scores <- GWAS(pheno_df, geno_df, plot = F, n.core = parallel::detectCores() - 2))
  
}



