library(aws.s3)
library(httr)
library(tidyverse)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
library(SNPRelate)
library(data.table)

# crop <- 'Corn'
# bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
# path <- '/H2H_Data/'

source('vaultCredentials.R')
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

geno <- s3read_using(fread, 
               object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/imputedGeno_wide.csv',
               bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

geno_all <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds',
                      bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

s3load(object = 'Corn/Corn_IDs.RData', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')

s3load(object = 'GPC_Data/Corn/Corn_GPC.RData',
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data')

CropIDs_pedigree <- CropIDs %>% 
  select(M.GERMPLASM.PEDIGREE, M.GERMPLASM.ORIGIN) %>% 
  unique()

geno <- geno[!duplicated(geno$Pedigree),]

snp_info <- s3read_using(fread, 
                        object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Archive/Sweetcorn20200703-1_MarkerMetadata.txt', 
                        bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

sample_id <- geno$Pedigree
snp_id <- colnames(geno)[-c(1:2)]
snp_pos <- snp_info$cM_times_1e6[order(match(snp_info$MRN, snp_id))]
snp_chr <- snp_info$LG[order(match(snp_info$MRN, snp_id))]
snp_allele <- paste(snp_info$allele0, snp_info$allele1, sep = '/')[order(match(snp_info$MRN, snp_id))]

# Create GDS file

snpgdsCreateGeno("/mnt/Corn/temperateCorn.gds", genmat = t(as.matrix(geno[,3:ncol(geno)])),
                 sample.id = sample_id, snp.id = snp_id,
                 snp.chromosome = snp_chr,
                 snp.position = snp_pos,
                 snp.allele = snp_allele, snpfirstdim=TRUE)

genofile <- snpgdsOpen('/mnt/Corn/temperateCorn.gds')

set.seed(1000)

snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2)
snpset.id <- unlist(unname(snpset))

# Run PCA
pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 4)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# 5.81 3.20 2.43 2.08 1.81 1.71 # first 5 principle component
sum(pc.percent, na.rm = T) # 34.714% variation

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

tab <- tab %>% 
  left_join(CropIDs, by = c('sample.id' = 'M.GERMPLASM.PEDIGREE'))

head(tab)

ggplot(tab, aes(x=EV2, y=EV1,  color=M.GERMPLASM.ORIGIN)) +
  geom_point() + 
  theme(legend.position = "none")


### Pick some of origins of interest

origin_2020pilot <- read.csv('/mnt/Corn/origin_2020pilot.csv')

snp_info <- s3read_using(fread, 
                         object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Archive/Sweetcorn20200703-1_MarkerMetadata.txt', 
                         bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
geno_all <- geno_all %>% 
  left_join(CropIDs_pedigree, by = c('Pedigree' = 'M.GERMPLASM.PEDIGREE'))

geno_sub <- geno_all %>% 
  filter(as.character(M.GERMPLASM.ORIGIN) %in% as.character(origin_2020pilot$origin))
length(unique(geno_sub$M.GERMPLASM.ORIGIN))

geno_sub <- geno_sub[!duplicated(geno_sub$Pedigree),]

sample_id <- geno_sub$Pedigree
snp_id <- colnames(geno_sub)[-c(1:2, ncol(geno_sub))]
snp_pos <- snp_info$cM_times_1e6[order(match(snp_info$MRN, snp_id))]
snp_chr <- snp_info$LG[order(match(snp_info$MRN, snp_id))]
snp_allele <- paste(snp_info$allele0, snp_info$allele1, sep = '/')[order(match(snp_info$MRN, snp_id))]


# Create GDS file

snpgdsCreateGeno("/mnt/Corn/temperateCorn_origin_2020pilot_1.gds", genmat = t(as.matrix(geno_sub[,3:(ncol(geno_sub) - 1)])),
                 sample.id = sample_id, snp.id = snp_id,
                 snp.chromosome = snp_chr,
                 snp.position = snp_pos,
                 snp.allele = snp_allele, snpfirstdim=TRUE)

genofile <- snpgdsOpen('/mnt/Corn/temperateCorn_origin_2020pilot_1.gds')

set.seed(1000)

snpset <- snpgdsLDpruning(genofile, ld.threshold = 0.2)
snpset.id <- unlist(unname(snpset))

# Run PCA
pca <- snpgdsPCA(genofile, snp.id = snpset.id, num.thread = 4)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# 6.06 4.05 3.43 2.82 2.67 2.15 # first 5 principle component
sum(pc.percent, na.rm = T) # 42.58% variation

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

tab <- tab %>% 
  left_join(CropIDs_pedigree, by = c('sample.id' = 'M.GERMPLASM.PEDIGREE')) %>% 
  left_join(GCP_TableView_noProg[, c('Pedigree', 'HETGP')], by = c('sample.id' = 'Pedigree'))

head(tab)

ggplot(tab, aes(x=EV2, y=EV1,  color=M.GERMPLASM.ORIGIN)) +
  geom_point() + 
  theme(legend.position = "none") + 
  ggtitle('Clustering by Origin: Origins in 2020 Pilot')

ggplot(tab, aes(x=EV2, y=EV1,  color=HETGP)) +
  geom_point() + 
  theme() + 
  ggtitle('Clustering by Heterotic: Origins in 2020 Pilot')

# geno_1 <- geno %>% 
#   select(-M.GERMPLASM.ORIGIN) %>% 
#   mutate(GermID = as.character(GermID)) %>% 
#   left_join(CropIDs_1, by = c('GermID' = 'M.GERMPLASM.X_ID'))
# 
# length(intersect(as.character(geno_1$M.GERMPLASM.ORIGIN), as.character(origin_2020pilot$origin)))
# intersect(as.character(geno_1$M.GERMPLASM.ORIGIN), as.character(origin_2020pilot$origin))
# # "SYW-6RLAC068/SHY-6S15-4568LN" "SHY-6SARH064/SHY-6RKCT030"    "SHY817-XT1130/SYW-6RLAC068"
# 
# 
# # Genotype from Implementation 
# 
# geno_implementation <- s3read_using(fread, 
#                                     object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/Geno_wide_SC1.csv',
#                                     bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
# 
# geno_all <- s3readRDS(object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/ImputedGeno_all.rds',
#                       bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')
#   
# length(intersect(as.character(geno_implementation$Pedigree), as.character(geno_all$Pedigree)))
# # Does geno_all include the implementation genotype? Yes, it does




# plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

# legend("bottomright", legend=levels(tab$M.GERMPLASM.ORIGIN), pch="o", col=1:nlevels(tab$M.GERMPLASM.ORIGIN))
