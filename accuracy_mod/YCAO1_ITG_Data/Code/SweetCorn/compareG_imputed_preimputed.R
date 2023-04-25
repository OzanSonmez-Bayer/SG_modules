# Compare non-imputed and imputed genotypic data

library(tidyverse)
library(dplyr)
library(aws.s3)
library(httr)

source('Credential/vaultCredentials.R')
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

origin_2020pilot <- read.csv('/mnt/Corn/origin_2020pilot.csv')

geno_all <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds',
                      bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

geno_preimputed <- s3readRDS(object = 'ycao1/ITG/Corn/gbs_preimpuatation_origin_2020.rds',
                             bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

geno_preimputed <- geno_preimputed[!duplicated(geno_preimputed$M.GERMPLASM.PEDIGREE),]


# Clean non-imputed data 

source('/repos/ITG_Data/Code/Shared/genoFunctions.R')

ids <- geno_preimputed %>% 
  select(M.GERMPLASM.PEDIGREE, sample)

geno_preimputed_1 <- geno_preimputed[, -c(1,2)]

geno_preimputed_1 <- apply(geno_preimputed_1, 2, as.character)
rownames(geno_preimputed_1) <- ids$M.GERMPLASM.PEDIGREE

geno_preimputed_2 <- genElement(geno_preimputed_1)
geno_preimputed_3 <- genClean(geno_preimputed_2)

rownames(geno_preimputed_3) <- ids$M.GERMPLASM.PEDIGREE

geno_preimputed_cln <- gen2Additive(geno_preimputed_3) # ACGT to 0,1,2
rownames(geno_preimputed_cln) <- ids$M.GERMPLASM.PEDIGREE


########### Compare Geno Before and After Imputation ################

overlapped_ped <- intersect(as.character(geno_all$Pedigree), rownames(geno_preimputed_cln))
overlapped_mrk <- intersect(colnames(geno_all), colnames(geno_preimputed_cln))

geno_imputed <- geno_all %>% 
  filter(as.character(Pedigree) %in% overlapped_ped)
geno_imputed <- geno_imputed[!duplicated(geno_imputed$Pedigree),]

rownames(geno_imputed) <- geno_imputed$Pedigree

geno_preimputed_cln <- geno_preimputed_cln[match(overlapped_ped, rownames(geno_preimputed_cln)), ]
geno_imputed <- geno_imputed[match(overlapped_ped, geno_imputed$Pedigree),]
geno_imputed_cln <- geno_imputed %>% 
  select(overlapped_mrk)

geno_preimputed_cln <- geno_preimputed_cln[, overlapped_mrk]

geno_preimputed_cln <- as.matrix(geno_preimputed_cln)
geno_imputed_cln <- as.matrix(geno_imputed_cln)

# switch 0 and 2
geno_preimputed_cln <- abs(geno_preimputed_cln - 2)

geno_compare <- rbind(geno_imputed_cln, geno_preimputed_cln)
geno_compare <- as.data.frame(geno_compare)

rownames(geno_compare) <- c(paste0(rownames(geno_preimputed_cln), '_imputed'),
                            paste0(rownames(geno_preimputed_cln), '_preimputed'))

geno_compare_t <- t(geno_compare)

# Calculate genetic distance 
g_distance <- gDist(geno_compare_t)

# Check up right matrix
g_distance_upright <- g_distance[1: nrow(g_distance)/2, 
                                 (ncol(g_distance)/2 + 1): ncol(g_distance)]
# check diag
range(diag(g_distance_upright))
hist(diag(g_distance_upright))

################### Compare Sister Lines within Origin 2020 Pilot ################

CropIDs_pedigree <- CropIDs %>% 
  select(M.GERMPLASM.PEDIGREE, M.GERMPLASM.ORIGIN) %>% 
  unique()

geno_all_1 <- geno_all[!duplicated(geno_all$Pedigree), ]

geno_all_1 <- geno_all_1 %>% 
  left_join(CropIDs_pedigree, by = c('Pedigree' = 'M.GERMPLASM.PEDIGREE'))

geno_2020 <- geno_all_1 %>% 
  filter(as.character(M.GERMPLASM.ORIGIN) %in% origin_2020pilot$origin)


geno_2020_1 <- geno_2020[,3:(ncol(geno_2020) - 1)]


geno_2020_1 <- t(geno_2020_1)
colnames(geno_2020_1) <- as.character(geno_2020$Pedigree)

sister_gDist <- gDist(geno_2020_1)

s3saveRDS(sister_gDist, object = 'ycao1/ITG/Corn/sister_gDist_2020.rds', 
          bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

dim(sister_gDist)

origin_2020pilot <- origin_2020pilot %>% 
  left_join(CropIDs_pedigree, by = c('origin' = 'M.GERMPLASM.ORIGIN'))

origin_2020pilot <- origin_2020pilot %>% 
  unique()

sister_gDist_sub <- sister_gDist[rownames(sister_gDist) %in% as.character(origin_2020pilot$M.GERMPLASM.PEDIGREE[which(origin_2020pilot$origin == 'SHW-6RAMQ040/SHY-6S18-3514GG')]),
                                 colnames(sister_gDist) %in% as.character(origin_2020pilot$M.GERMPLASM.PEDIGREE[which(origin_2020pilot$origin == 'SHW-6RAMQ040/SHY-6S18-3514GG')])]

range(sister_gDist_sub, na.rm = T) # 0.03900022 0.20079724
# install snpReady

# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("impute")
# 
# install.packages('snpReady')
# 
# library(snpReady)
# 
# geno_preimputed_cln <- raw.data(data = as.matrix(geno_preimputed_3), 
#                    frame = 'wide',
#                    base = TRUE, 
#                    sweep.sample = 0.5,
#                    call.rate = 0.95,
#                    maf = 0.05, 
#                    imput = TRUE, 
#                    imput.type = 'wright',
#                    outfile = '012')

#################### Pairwise Comparison ####################
geno_preimputed_012 <- geno_preimputed_2

for (i in 1:ncol(geno_preimputed_2)){
  for (j in 1:nrow(geno_preimputed_2)){
    if (as.character(geno_preimputed_2[j, i]) %in% c('AA', 'CC')){
      geno_preimputed_012[j,i] <- 2
    }else if (as.character(geno_preimputed_2[j, i]) %in% c('GG', 'TT')){
      geno_preimputed_012[j,i] <- 0
    }else if (as.character(geno_preimputed_2[j, i]) == 'NA'){
      geno_preimputed_012[j,i] <- NA
    }else if (!as.character(geno_preimputed_2[j, i]) %in% c('AA', 'CC', 'GG', 'TT', 'NA'))
      geno_preimputed_012[j,i] <- 1
  }
}

geno_preimputed_012 <- apply(geno_preimputed_012, MARGIN = 2, FUN = as.numeric)

rownames(geno_preimputed_012) <- rownames(geno_preimputed_2)

geno_preimputed_012[1:5,1:5]

# Pairwise copmarision

dim(geno_compare_t)

perct_diff_true <- matrix(data = NA, nrow = 1444, ncol = 3)


for (tt in 1:(ncol(geno_compare_t)/2)){
  temp <- unname(table(geno_compare_t[, tt] == geno_compare_t[, (tt + ncol(geno_compare_t)/2)]))/901
  perct_diff_true[tt, ] <- c(temp, 1-sum(temp))
}

perct_df <- data.frame(Pedigree = colnames(geno_compare_t)[1:1444], 
                       Perct_Diff = perct_diff_true[,1], 
                       Perct_Same = perct_diff_true[,2],
                       Perct_Miss = perct_diff_true[,3])

perct_df$Pedigree <- str_remove(perct_df$Pedigree, '_imputed')

s3write_using(perct_df, FUN = write.csv, object = 'ycao1/ITG/Corn/QAQC_Geno/check_origin_2020.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


############# Phenotypic Distance #################

library(data.table)

pheno <- s3readRDS(object = 'ycao1/ITG/Corn/pheno_origin_hbc.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

# SSGBLUP-GCA data 

gca <- s3read_using(fread, object = 'ycao1/ITG/Corn/SSGBLUP_GCA_masking_2020_02222021.csv',
                    bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

gca_2020 <- gca %>% 
  filter(Origin %in% origin_2020pilot$origin) # only interested in the lines from 2020 pilot

gca_2020_mat <- as.matrix(gca_2020[, -c(1:4)])

# gca_2020_mat <- apply(gca_2020_mat, MARGIN = 2, FUN = function(x){x[which(is.na(x))] = mean(x, na.rm = T)})

for (i in 1:ncol(gca_2020_mat)){
  gca_2020_mat[,i][which(is.na(gca_2020_mat[,i]))] <- mean(gca_2020_mat[,i], na.rm = T)
  gca_2020_mat[,i] <- (gca_2020_mat[,i] - mean(gca_2020_mat[,i]))/sqrt(var(gca_2020_mat[,i]))
}

# gca_2020_mat[,1][which(is.na(gca_2020_mat[,1]))] <- mean(gca_2020_mat[,1], na.rm = T)
# gca_2020_mat[,2][which(is.na(gca_2020_mat[,2]))] <- mean(gca_2020_mat[,2], na.rm = T)
# gca_2020_mat[,3][which(is.na(gca_2020_mat[,3]))] <- mean(gca_2020_mat[,3], na.rm = T)
# 
# gca_2020_mat <- as.data.frame(gca_2020_mat)

pheno_dist <- dist(gca_2020_mat, method = 'euclidean')
pheno_dist <- as.matrix(pheno_dist)
colnames(pheno_dist) <- as.character(gca_2020$PEDIGREE_NAME)
rownames(pheno_dist) <- as.character(gca_2020$PEDIGREE_NAME)

gDist_2020 <- s3readRDS(object = 'ycao1/ITG/Corn/sister_gDist_2020.rds', 
                        bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

# Calculate Euclidean distance with genotypic data

gDist_e <- dist(geno_2020_1, method = 'euclidean')
gDist_e <- as.matrix(gDist_e)

# pick one origin 

cor_df <- data.frame(origin = unique(gca_2020$Origin), correlation = NA)

for(i in unique(gca_2020$Origin)){
  ped <- as.character(gca_2020$PEDIGREE_NAME[gca_2020$Origin == i])
  
  pDist_sub <- pheno_dist[ped, ped]
  gDist_sub <- gDist_2020[rownames(gDist_2020) %in% ped, colnames(gDist_2020) %in% ped]
  
  
  pDist_sub <- pDist_sub[rownames(gDist_sub), colnames(gDist_sub)]
  pDist_sub <- pDist_sub[match(rownames(gDist_sub), rownames(pDist_sub)), match(rownames(gDist_sub), rownames(pDist_sub))]
  
  test <- pDist_sub[upper.tri(pDist_sub)]
  test_1 <- gDist_sub[upper.tri(gDist_sub)]
  
  cor(test, test_1) # low correlation, could be different calculation methods
  
  # what if we normalize it
  test_norm <- (test - mean(test))/sqrt(var(test))
  test_1_norm <- (test_1 - mean(test_1))/sqrt(var(test_1))
  
  temp <- cor(test_norm, test_1_norm) # No, still the same. 
  
  cor_df[cor_df$origin == i, 'correlation'] <- temp
}


ggplot(data = cor_df,aes(y=correlation, x=origin)) + 
  geom_point(color ='red') + 
  geom_text_repel(aes(label = cor_df$origin), size = 3.5) + 
  theme_classic(base_size = 10) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
