# SC1 population

library(readxl)
library(aws.s3)
library(httr)
library(dplyr)
library(tidyverse)

sc1_sets <- xlsx::read.xlsx('Corn/SC1_info.xlsx', sheetIndex = 1)
sc1_fields <- xlsx::read.xlsx('Corn/SC1_info.xlsx', sheetIndex = 2)

source('vaultCredentials.R')
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

path <- '/H2H_Data/'
crop <- 'Corn'
infile <- paste0(path,crop,"/",'2020',"/AllYearData.RData")
s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")

SC1_data <- YearPhenos %>% 
  filter(GROWSEASON == '2020:04'&
         TEST_SET_NAME %in% as.character(sc1_sets[,1]) &
         FIELD_NAME %in% as.character(sc1_fields[,1])) 

pheno_d=subset(SC1_data,DESCRIPTOR %in% c('PUCCSO', 'SETOTU'))
SC1_data[rownames(pheno_d),]$OBSRVTN_REF_CD=paste0(pheno_d$OBSRVTN_REF_CD,'_',pheno_d$DESCRIPTOR)

SC1_data <- SC1_data %>% 
  filter(OBSRVTN_REF_CD %in% c('EDIA', 'ELEN', 'ESHK', 'HSC', 'PA', 'QUAL', 'RTLR', 'S50', 
                               'S50_BE', 'SC_FG', 'SC_FR', 'SC_HE', 'SC_PR', 'SC_SE', 'SC_TF', 'SDV',
                               'TILLR', 'LDSR_PUCCSO', 'LDSR_SETOTU'))

pheno_d_1 <- subset(itg_data,DESCRIPTOR %in% c('PUCCSO', 'SETOTU'))
itg_data[rownames(pheno_d_1),]$OBSRVTN_REF_CD=paste0(pheno_d_1$OBSRVTN_REF_CD,'_',pheno_d_1$DESCRIPTOR)

pheno <- rbind(itg_data, SC1_data)
pheno$TREP <- paste(pheno$REP_NUMBER,pheno$TRACK_ID,sep='_')

s3saveRDS(pheno, object = '/ycao1/ITG/pheno.rds', 
          bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

# 
# %>%
#   mutate(P2 = if_else(is.na(P2), P1, P2)) %>%
#   select(PEDIGREE_NAME, P1, P2) %>%
#   distinct()


GT <- readxl::read_excel('Corn/SweetCorn_SC1_GenotypingStatus.xlsx')
GT_unique <- GT %>% 
  select(HybridPedigree, FatherPedigree, MotherPedigree) %>% 
  distinct()

# Comparison 

hybrid_h2h <- unique(SC1_data$PEDIGREE_NAME)
hybrid_genomic <- unique(GT_unique$HybridPedigree)

length(hybrid_h2h)
length(hybrid_genomic)

length(intersect(hybrid_h2h, hybrid_genomic))

# parents 

parent_h2h <- unique(c(SC1_data$P1, SC1_data$P2))
parent_genomic <- unique(c(GT_unique$FatherPedigree, GT_unique$MotherPedigree))

length(parent_h2h)
length(parent_genomic)

length(intersect(parent_h2h, parent_genomic))
