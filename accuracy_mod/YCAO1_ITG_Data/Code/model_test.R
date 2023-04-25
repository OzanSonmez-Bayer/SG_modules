library(asreml)
library(asremlPlus)
library(tidyverse)
library(aws.s3)
library(httr)

install.packages('/mnt/packages/parallelly_1.20.0.tar.gz', repos = NULL, type = 'source')
library(parallelly)

install.packages('/mnt/packages/doBy_4.6-2.tar.gz', repos = NULL, type = 'source')
library(doBy)

source('Credential/vaultCredentials.R')
source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

pheno <- s3readRDS(object = 'ycao1/ITG/Corn/pheno.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


dat <- pheno %>% 
  filter(OBSRVTN_REF_CD == 'ELEN') %>% 
  filter(GROWSEASON %in% c("2020:04")) %>% 
  mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE))

# dat_origin <- dat %>% 
  # filter(grepl('6S_SEW-6S17-4402ZZ/SHY-6S18-0773EE', P1))

dat_origin <- dat


############## H2H BLUE ##########################

# blue <- blue_fun(dataset = dat, crop = 'Corn')

if(length(grep("Windows",Sys.getenv("OS")))>0) {
  NumberOfCores = 1
} else {
  NumberOfCores = availableCores()/2
  
  print(paste0("Number of cores detected ",availableCores()))
  
}

if(NumberOfCores < 1) { NumberOfCores = 1}


NumberOfCrossNames <- length(unique(dat$CROSS_NAME))

CrossesPerRun <- ceiling(NumberOfCrossNames/10)

if(CrossesPerRun > 1500) { 
  CrossesPerRun = 1500 
}else if (CrossesPerRun < 500) { 
    CrossesPerRun = NumberOfCrossNames 
    }

ListCrossNames <- unique(dat$CROSS_NAME)
ListOfList <- split(ListCrossNames,as.numeric(gl(length(ListCrossNames),CrossesPerRun,length(ListCrossNames))))

ListOf_PredBlue <- parallel::mclapply(ListOfList, function(x) blue_fun(dat,x), mc.cores = NumberOfCores, mc.silent = FALSE)
pred_blu <- do.call(rbind,ListOf_PredBlue)

############### H2H PBLUP ########################

mdl_h2h <- pblup_fun(dataset = dat, crop = 'Corn') 

sca_h2h <- mdl_h2h[[2]]
gca_h2h <- mdl_h2h[[1]]

############### New BLUP Model #####################

mdl_new_gca <- asreml(fixed = TRAIT_VALUE ~ 1 + GROWSEASON + REPETITION + BR_FIELD_ID:TREP, 
                  random = ~ P1 + P2, 
                  data = dat, 
                  na.method.X = 'include', 
                  na.method.Y = 'include')

pred_new <- predict(mdl_new_gca, classify = 'P1')
P1_gca <- pred_new$predictions$pvals
colnames(P1_gca)[1] <- 'Parent'

head(P1_gca)

P2_gca <- predict(mdl_new_gca, classify = 'P2')$predictions$pvals
colnames(P2_gca)[1] <- 'Parent'

head(P2_gca)

gca_new <- rbind(P1_gca, P2_gca) # some parents are being used as female and male

gca_new[duplicated(gca_new$Parent),]

gca_new <- gca_new %>% 
  group_by(Parent) %>% 
  summarise(predicted_value = mean(predicted.value)) %>% 
  as.data.frame()

head(gca_new)


######################### GCA model with GINV #############################

Crop <- 'Corn'
source('/mnt/Credential/Credential_sdi.R')

Sys.setenv("CLIENT_ID" = CLIENT_ID,
           "CLIENT_SECRET" = CLIENT_SECRET)

source('/repos/ITG_Data/Code/Shared/deepped.R')



ped_deep <- dat %>% 
  select(PEDIGREE_NAME, P1, P2, LTYPE) %>% 
  dplyr::rename(pedigree = PEDIGREE_NAME) %>% 
  as.data.frame()

ped_deep <- ped_deep[!duplicated(ped_deep$pedigree),]

head(ped_deep)

female_parent <- unique(ped_deep$P1)
male_parent <- unique(ped_deep$P2)

female_germid <- CropIDs %>% 
  dplyr::select(M.GERMPLASM.X_ID, M.GERMPLASM.PEDIGREE) %>% 
  filter(M.GERMPLASM.PEDIGREE %in% female_parent)

male_germid <- CropIDs %>% 
  dplyr::select(M.GERMPLASM.X_ID, M.GERMPLASM.PEDIGREE) %>% 
  filter(M.GERMPLASM.PEDIGREE %in% male_parent)

ped_file <- make_ped_file(ped_deep, crop_ids = CropIDs, n_gen = 3)

female_ped <- ped_file$ped %>% 
  filter(ID %in% female_germid$M.GERMPLASM.X_ID) %>% 
  left_join(ped_file$distinct_lines, by = 'ID')

female_ped <- female_ped[!duplicated(female_ped$ID),]

male_ped <- ped_file$ped %>% 
  filter(ID %in% male_germid$M.GERMPLASM.X_ID)


################# One Origin ###########

ped_deep_sub <- dat_origin %>% 
  select(PEDIGREE_NAME, P1, P2, LTYPE) %>% 
  dplyr::rename(pedigree = PEDIGREE_NAME) %>% 
  as.data.frame()


ped_deep_1 <- ped_deep_sub[!duplicated(ped_deep$pedigree),]

head(ped_deep_1)

female_parent <- unique(ped_deep_1$P1)
male_parent <- unique(ped_deep_1$P2)

female_germid <- CropIDs %>% 
  dplyr::select(M.GERMPLASM.X_ID, M.GERMPLASM.PEDIGREE) %>% 
  filter(M.GERMPLASM.PEDIGREE %in% female_parent) %>% 
  unique()

male_germid <- CropIDs %>% 
  dplyr::select(M.GERMPLASM.X_ID, M.GERMPLASM.PEDIGREE) %>% 
  filter(M.GERMPLASM.PEDIGREE %in% male_parent) %>% 
  unique()

ped_file_1 <- make_ped_file(ped_deep_sub, crop_ids = CropIDs, n_gen = 3)

ginv_all <- asreml.Ainverse(ped_file_1$ped, fgen = list("inbreeding", 5))



female_ped <- ped_file_1$ped %>%
  filter(ID %in% female_germid$M.GERMPLASM.X_ID) 

female_ped <- female_ped[!duplicated(female_ped$ID),]

male_ped <- ped_file_1$ped %>%
  filter(ID %in% male_germid$M.GERMPLASM.X_ID)


female_ginv <- asreml.Ainverse(female_ped, fgen = list("inbreeding", 5))
male_ginv <- asreml.Ainverse(male_ped, fgen = list('inbreeding', 5))

dat_origin_1 <- dat_origin %>% 
  left_join(ped_file_1$distinct_lines[ , c("ID", "germplasm_id", "P1")], by = "P1") %>% 
  dplyr::rename(P1_ID = ID) %>% 
  left_join(ped_file_1$distinct_lines[ , c("ID", "germplasm_id", "P2")], by = 'P2') %>% 
  dplyr::rename(P2_ID = ID) %>% 
  select(TRAIT_VALUE, P1_ID, P2_ID, GROWSEASON, FIELD_NAME, TREP) 

dat_origin_1 <- dat_origin_1 %>% 
  unique()

dim(dat_origin_1)
head(dat_origin_1)


mdl_new_gca_ginv <- asreml(fixed = TRAIT_VALUE ~ 1, 
                      random = ~ P1_ID + P2_ID, 
                      data = dat_origin_1, 
                      na.method.X = 'include', 
                      na.method.Y = 'include',
                      ginverse = list(P1_ID = female_ginv$ginv, P2_ID = male_ginv$ginv),
                      workspace = 320e+6, pworkspace = 320e+6)

# Singularity in the A mat?

female_A <- asreml.sparse2mat(female_ginv$ginv)
female_A <- round(solve(female_A),2) 

male_A <- asreml.sparse2mat(male_ginv$ginv)
male_A <- round(solve(male_A),2) 

pred_new <- predict(mdl_new_gca_ginv, classify = 'P1_ID')




################# Form Ginverse with genotypic data 

geno <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds', 
                  bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

female_geno <- geno %>% 
  filter(Pedigree %in% female_parent)

female_geno_mat <- female_geno[, 3:ncol(female_geno)] %>% 
  as.matrix()

rownames(female_geno_mat) <- female_geno$Pedigree


male_geno <- geno %>% 
  filter(Pedigree %in% male_parent)

male_geno_mat <- male_geno[, 3:ncol(male_geno)] %>% 
  as.matrix()

rownames(male_geno_mat) <- male_geno$Pedigree


install.packages('/mnt/packages/AGHmatrix_0.0.4.tar.gz', repos = NULL, type = 'source')
library(AGHmatrix)

male_gmat <- Gmatrix(male_geno_mat, method = 'VanRaden')
male_gmat_inv <- solve(male_gmat, tol = 1e-17)

male_gmat_asreml <- formatmatrix(male_gmat_inv,return = TRUE)
colnames(male_gmat_asreml) <- c('row', 'column', 'Ainverse')
rownames(male_gmat_asreml) <- 1:nrow(male_gmat_asreml)

head(male_gmat_asreml)


female_gmat <- Gmatrix(female_geno_mat, method = 'VanRaden')
female_gmat_inv <- solve(female_gmat, tol = 1e-20)
female_gmat_asreml <- formatmatrix(female_gmat_inv, return = TRUE)
colnames(female_gmat_asreml) <- c('row', 'column', 'Ainverse')
rownames(female_gmat_asreml) <- 1:nrow(female_gmat_asreml)


# library(snpReady)
# 
# male_gmat <- G.matrix(male_geno_mat, method = 'VanRaden', format = 'long')

mdl_new_gca_ginv <- asreml(fixed = TRAIT_VALUE ~ 1, 
                           random = ~ P1 + P2, 
                           data = dat_origin, 
                           na.method.X = 'include', 
                           na.method.Y = 'include',
                           ginverse = list(P1 = female_gmat_asreml , P2 = male_gmat_asreml),
                           workspace = 320e+6, pworkspace = 320e+6)


########################## Compare GCA between Two Models #################

gca_compare <- gca_h2h %>% 
  left_join(gca_new, by = c('PEDIGREE_NAME' = 'Parent'))

plot(gca_compare$predicted.value, gca_compare$predicted_value)

cor(gca_compare$predicted.value, gca_compare$predicted_value, method = 'spearman')

# Pearson = 0.90, Spearman = 0.89


############ New BLUP Model for SCA #################

mdl_new_sca <- asreml(fixed = TRAIT_VALUE ~ 1 + GROWSEASON + REPETITION + BR_FIELD_ID:TREP, 
                      random = ~ PEDIGREE_NAME, 
                      data = dat, 
                      na.method.X = 'include')

pred_sca <- predict(mdl_new_sca, classify = 'PEDIGREE_NAME')$predictions$pvals

head(pred_sca)


########### Compare SCA between Two Models ###############

sca_compare <- sca_h2h %>% 
  inner_join(pred_sca, by = 'PEDIGREE_NAME')

plot(sca_compare$predicted.value.x, sca_compare$predicted.value.y)

cor(sca_compare$predicted.value.x, sca_compare$predicted.value.y, method = 'spearman')

# spearmean = 0.83, pearson = 0.84


# raw data distribution 

hist(as.numeric(dat$TRAIT_VALUE))
hist(pred_sca$predicted.value)


########## Top 5% ############

top10 <- sca_h2h %>% 
  arrange(desc(predicted.value)) %>% 
  # slice_head(prop = 0.05)
  slice_head(n=10) %>% 
  select(PEDIGREE_NAME)

top10_1 <- pred_sca %>% 
  arrange(desc(predicted.value)) %>% 
  slice_head(n=10) %>% 
  select(PEDIGREE_NAME)

top10
top10_1

intersect(top10, top10_1)
