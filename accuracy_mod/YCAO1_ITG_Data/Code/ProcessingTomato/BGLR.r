# BGLR
library(aws.s3)
library(httr)
install.packages('BGLR')
library(BGLR)
library(tidyverse)

source('vaultCredentials.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

geno <- s3readRDS(object = 'Tomato/ProcessingTomato/9Z/geno_all.rds', 
                  bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

pheno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/pheno_training.rds',
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

### Center G matrix 

geno_centered <- geno[, 3:ncol(geno)]
rownames(geno_centered) <- geno$Pedigree
geno_centered <- geno_centered - 1
geno_centered <- scale(geno_centered, center = TRUE, scale = TRUE)
G <- tcrossprod(geno_centered)/ncol(geno_centered)

# Computing the eigen-value
EVD <- eigen(G)
saveRDS(EVD, file = '/mnt/GWS/ProcTomato/BGLR/EVD.rds')

### Create a A matrix in Asreml 

library(asreml)
library(asremlPlus)


Crop <- 'Tomato'
trait <- 'AVGHBRIX'
n_gen <- 5

source('/mnt/Credential_sdi.R')
source('/repos/ITG_Data/Code/Shared/deepped.R')


pheno_sub <- pheno %>% 
  filter(OBSRVTN_REF_CD == trait) %>% 
  rename(pedigree = PEDIGREE_NAME) %>% 
  mutate(TREP = paste(REP_NUMBER, TRACK_ID,sep='_'),
         TRAIT_VALUE = as.numeric(TRAIT_VALUE))

# mask testing season
yNA <- pheno_sub$TRAIT_VALUE
yNA[which(pheno_sub$GROWSEASON == '2020:04')] <- NA


# pd=unique(data.frame(bdat$CROSS_NAME,bdat$P1,bdat$P2))
# pd=pd[!duplicated(pd[,1]),]
# names(pd)=c('CROSS_NAME','P1','P2')
# pd$Selfing=rep(0,rep=nrow(pd))
# pd.ainv=asreml.Ainverse(pd,fgen = c('Selfing',5))

new_token(fingerprint = F)
gen5 <- make_ped_file(pheno_sub, crop_ids = CropIDs, n_gen = n_gen)
A_matrix5 <- asreml.Ainverse(gen5$ped, fgen = list("inbreeding", n_gen))
Amat <- round(solve(asreml.sparse2mat(A_matrix5$ginv)),2) 

colnames(Amat) <- A_matrix5$pedigree$PEDIGREE_NAME
rownames(Amat) <- A_matrix5$pedigree$PEDIGREE_NAME


ETA <- list(list(~factor(GROWSEASON) + factor(FIELD_NAME) + factor(TREP), data = pheno_sub, model = 'BRR'),
            list(K = Amat, model = 'RKHS'),
            list(V = EVD$vectors, d = EVD$values, model = 'RKHS'))

fm <- BGLR(y = yNA, ETA = ETA, nIter = 12000, burnIn = 2000, saveAt = '/mnt/GWS/Corn/BGLR/' )
