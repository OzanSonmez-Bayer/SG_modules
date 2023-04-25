################################# GBLUP with Lower Marker Density #######################

################ LD prunning 

# library(devtools)
# install_github("zhengxwen/gdsfmt")
# install_github("zhengxwen/SNPRelate")

library(gdsfmt)
library(SNPRelate)
library(jsonlite)
library(data.table)
library(Rcpp)
library(aws.s3)
library(tidyverse)
library(rrBLUP)
library(httr)

source('vaultCredentials.R')

bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
crop <- 'Tomato'

# Credentials

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

geno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/geno_all.rds', 
                  bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

# Import marker metadata

marker_meta <- s3read_using(fread, object = 'Tomato/ProcessingTomato/9Z/ProcessingTomato20200615-1_MarkerMetadata.txt', 
                            bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

geno <- geno[!duplicated(geno$Pedigree),]
sample_id <- geno$Pedigree
snp_id <- colnames(geno)[-c(1:2)]
snp_pos <- marker_meta$cM_times_1e6[order(match(marker_meta$MRN, snp_id))]
snp_chr <- marker_meta$LG[order(match(marker_meta$MRN, snp_id))]
snp_allele <- paste(marker_meta$allele0, marker_meta$allele1, sep = '/')[order(match(marker_meta$MRN, snp_id))]

# Create GDS file

# snpgdsCreateGeno("ProcessingTomato.gds", genmat = t(as.matrix(geno[,3:ncol(geno)])),
#                  sample.id = sample_id, snp.id = snp_id,
#                  snp.chromosome = snp_chr,
#                  snp.position = snp_pos,
#                  snp.allele = snp_allele, snpfirstdim=TRUE)
# 
# genofile <- snpgdsOpen("ProcessingTomato.gds")


set.seed(20200804)

# snpset <- snpgdsLDpruning(genofile, maf = 0.05, ld.threshold=0.05)
# 
# geno_prunned <- geno %>% 
#   select(unname(unlist(snpset)))


######### GBLUP with rrBLUP ################



pheno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/pheno_training.rds',
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

trait_list <- c('AVGHBRIX', 'AVJB', 'FQUAL', 'FRFRMH', 'FZUNI', 'LBRIX', 'MAT', 'MATR', 'OST', 'PLCMP', 'PLTVG', 'PYLDPA', 'SAMPUW')

accuracy <- data.frame(trait = NA, training = NA, testing = NA)

for (trait in trait_list){
  
  print(trait)
  
  pheno_sub <- pheno %>% 
    filter(OBSRVTN_REF_CD == trait) %>% 
    mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE))
  
  # individuals from testing population 
  ped_testing <- pheno_sub %>% 
    filter(GROWSEASON == '2019:04') %>% 
    select(PEDIGREE_NAME) %>% 
    distinct()
  
  ped_testing <- ped_testing[,'PEDIGREE_NAME']
  
  # LSmean
  
  pheno_dat <- pheno_sub %>% 
    group_by(PEDIGREE_NAME) %>% 
    summarise(lsmean = mean(TRAIT_VALUE, na.rm = T)) %>% 
    mutate(lsmean_observed = lsmean)
  
  # Mask phenotypic value from testing population 
  pheno_dat$lsmean[pheno_dat$PEDIGREE_NAME  %in% ped_testing] <- NA
  
  ### Match pheno with genotype
  
  # geno_prunned$Pedigree <- geno$Pedigree
  
  
  # pheno_geno <- pheno_dat %>%
  #   inner_join(geno_prunned, by = c('PEDIGREE_NAME' = 'Pedigree'))
  
  # pheno_geno <- pheno_dat %>%
  #   inner_join(geno, by = c('PEDIGREE_NAME' = 'Pedigree'))
  # 
  pheno_geno <- pheno_dat %>%
    full_join(geno, by = c('PEDIGREE_NAME' = 'Pedigree'))
  
  geno_cln <- pheno_geno[,colnames(pheno_geno)[which(!colnames(pheno_geno) %in% c('PEDIGREE_NAME', 'lsmean', 'lsmean_observed', 'GermID'))]] - 1
  
  #### Kinship Matrix
  kinship_mat <- A.mat(as.matrix(geno_cln), return.imputed = T)  
  
  #### GBLUP 
  rrblup_mdl <- mixed.solve(y = pheno_geno$lsmean, K = kinship_mat$A)
  
  gblup <- as.vector(rrblup_mdl$u) + as.vector(rrblup_mdl$beta)
  # 
  cor_train <- cor(pheno_geno$lsmean_observed[which(!is.na(pheno_geno$lsmean))], gblup[which(!is.na(pheno_geno$lsmean))])
  # 
  cor_test <- cor(pheno_geno$lsmean_observed[which(is.na(pheno_geno$lsmean))], gblup[which(is.na(pheno_geno$lsmean))])
  
  accuracy <- rbind(accuracy, data.frame(trait = trait, training = cor_train, testing = cor_test))
  
  saveRDS(gblup, file = paste0('ProcTomato/rrBLUP/', trait, '.RDS'))
  
}

################## Test Out different LD Thresholds ################



trait_list <- c('AVGHBRIX', 'AVJB', 'FQUAL', 'FRFRMH', 'FZUNI', 'LBRIX', 'MAT', 'MATR', 'OST', 'PLCMP', 'PLTVG', 'PYLDPA', 'SAMPUW')


accuracy <- data.frame(trait = NA, marker_num = NA, training = NA, testing = NA)

for (trait in trait_list){
  
  print(trait)
  # trait <- 'ELEN'
  
  ld_threshold_ls <- c(0.05, 0.2, 0.5, 0.7, 0.9)
  
  # accuracy <- data.frame(trait = trait, marker_num = NA, training = NA, testing = NA)
  
  for (ld in ld_threshold_ls){
    
    print(ld)
    
    set.seed(20200804)
    
    snpset <- snpgdsLDpruning(genofile, maf = 0.05, ld.threshold = ld)
    
    geno_prunned <- geno %>% 
      select(unname(unlist(snpset)))
    
    
    pheno_sub <- pheno %>% 
      filter(OBSRVTN_REF_CD == trait) %>% 
      mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE))
    
    # individuals from testing population 
    ped_testing <- pheno_sub %>% 
      filter(GROWSEASON == '2019:04') %>% 
      select(PEDIGREE_NAME) %>% 
      distinct()
    
    ped_testing <- ped_testing[,'PEDIGREE_NAME']
    
    # LSmean
    
    pheno_dat <- pheno_sub %>% 
      group_by(PEDIGREE_NAME) %>% 
      summarise(lsmean = mean(TRAIT_VALUE, na.rm = T)) %>% 
      mutate(lsmean_observed = lsmean)
    
    # Mask phenotypic value from testing population 
    pheno_dat$lsmean[pheno_dat$PEDIGREE_NAME  %in% ped_testing] <- NA
    
    ### Match pheno with genotype
    
    geno_prunned$Pedigree <- geno$Pedigree
    
    
    pheno_geno <- pheno_dat %>%
      inner_join(geno_prunned, by = c('PEDIGREE_NAME' = 'Pedigree'))
    
    geno_cln <- pheno_geno[,colnames(pheno_geno)[which(!colnames(pheno_geno) %in% c('PEDIGREE_NAME', 'lsmean', 'lsmean_observed', 'GermID'))]] - 1
    
    #### Kinship Matrix
    kinship_mat <- A.mat(as.matrix(geno_cln), return.imputed = F, n.core = 60)  
    
    #### GBLUP 
    rrblup_mdl <- mixed.solve(y = pheno_geno$lsmean, K = kinship_mat)
    
    gblup <- as.vector(rrblup_mdl$u) + as.vector(rrblup_mdl$beta)
    # 
    cor_train <- cor(pheno_geno$lsmean_observed[which(!is.na(pheno_geno$lsmean))], gblup[which(!is.na(pheno_geno$lsmean))])
    # 
    cor_test <- cor(pheno_geno$lsmean_observed[which(is.na(pheno_geno$lsmean))], gblup[which(is.na(pheno_geno$lsmean))])
    
    accuracy <- rbind(accuracy, data.frame(trait = trait, marker_num = ncol(geno_prunned), training = cor_train, testing = cor_test))
    
    saveRDS(gblup, file = paste0('ProcTomato/rrBLUP/', trait, '_', ld, '.RDS'))
    
  }
  
  # accuracy <- accuracy[-1,]
  # accuracy$trait <- trait
  # accuracy$LD.threshold <- ld_threshold_ls
  
}


###################### Get GBLUP from Asreml ###################

# Kinship matrix including parental lines 

dim(geno)

G <- A.mat(as.matrix(geno[, 3:ncol(geno)]), return.imputed = T,n.core = 60)
G_mat <- G$A

colnames(G_mat) <- geno$Pedigree
rownames(G_mat) <- geno$Pedigree

realizedPD <- nearPD(G_mat, keepDiag = T)
G_base <- matrix(realizedPD[[1]]@x, nrow = realizedPD[[1]]@Dim[1])
G_base <- G_base + diag(0.01, nrow = nrow(G_base))
attr(G_base, 'dimnames') <- realizedPD[[1]]@Dimnames
class(G_base) <- 'relationshipMatrix'
str(G_base)
summary(G_base)

G_inv <- write.relationshipMatrix(G_base, sorting = 'ASReml', type = 'ginv')
head(attr(G_inv, 'rowNames'))
names(G_inv) <- c('row', 'column', 'coefficient')
head(G_inv)


# Mask phenotypic data from 2020:04

pheno$TRAIT_VALUE[which(pheno$GROWSEASON == '2020:04')] <- NA

GBLUP_asreml <- asreml(fixed = TRAIT_VALUE ~ GROWSEASON + BR_FIELD_ID:TREP + REPEITTION,
                       random = ~giv(PEDIGREE_NAME, var = T), 
                       ginverse = list(PEDIGREE_NAME = G_inv), 
                       rcov = ~units, 
                       data = pheno, 
                       na.method.X = 'include',
                       verbose = TRUE,
                       workspace = 320e+6, pworkspace = 320e+6)