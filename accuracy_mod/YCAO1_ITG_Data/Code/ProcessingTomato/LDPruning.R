# LD Prunning

library(devtools)
install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")

library(gdsfmt)
library(SNPRelate)
library(jsonlite)
library(data.table)
library(Rcpp)
library(aws.s3)
library(tidyverse)

source('Code/Shared/misFunctions.R')

AWS.Authorization('ycao1')

path <- 'ycao1/ITG_data/ProcTomato/'
bucketname <- 'genome-analytics-perm-space'
geno <- s3read_using(fread, 
                     object = 's3://genome-analytics-perm-space/ycao1/ITG_data/ProcTomato/Geno/imputedGeno.csv')


snp_info <- s3read_using(fread, 
                         object = 's3://genome-analytics-perm-space/shared/Veg_Genotypes/GWS_Pilots/Tomato/ProcessingTomato/ProcessingTomato20200615-1_MarkerMetadata.txt')


sample_id <- geno$Pedigree
snp_id <- colnames(geno)[-c(1:2)]
snp_pos <- snp_info$cM_times_1e6[order(match(snp_info$MRN, snp_id))]
snp_chr <- snp_info$LG[order(match(snp_info$MRN, snp_id))]
snp_allele <- paste(snp_info$allele0, snp_info$allele1, sep = '/')[order(match(snp_info$MRN, snp_id))]

# Create GDS file

snpgdsCreateGeno("Proctomato.gds", genmat = t(as.matrix(geno[,3:ncol(geno)])),
                 sample.id = sample_id, snp.id = snp_id,
                 snp.chromosome = snp_chr,
                 snp.position = snp_pos,
                 snp.allele = snp_allele, snpfirstdim=TRUE)

genofile <- snpgdsOpen("Proctomato.gds")

set.seed(20200804)

snpset <- snpgdsLDpruning(genofile, maf = 0.05, ld.threshold=0.2)

geno_prunned <- geno %>% 
  select(unname(unlist(snpset)))

