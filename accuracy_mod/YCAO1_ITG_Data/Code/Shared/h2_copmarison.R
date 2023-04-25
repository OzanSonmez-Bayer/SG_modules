# Checking heritability

library(aws.s3)
library(tidyverse)
library(ggplot2)
library(ggrepel)

source('credential.R')
source('misFunctions.R')

############################## Tomato ##################################
path <- ('ycao1/ITG_data/ProcTomato/')
bucketname <- 'genome-analytics-perm-space'

dfBucket_blup <- get_bucket_df(bucketname, paste0(path, 'BLUP'))
file_list_blup <- dfBucket_blup$Key[grep('h2', dfBucket_blup$Key)]

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/h2_PCM0.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/h2_PCM0_deep.csv"
heritabilityVis(file1, file2, "PCM0: Heritability Comparison")

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/h2_PCM1.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/h2_PCM1_deep.csv"
heritabilityVis(file1, file2)


########################### Sweet Corn #################################
path <- ('ycao1/ITG_data/SweetCorn/')
bucketname <- 'genome-analytics-perm-space'

dfBucket_blup <- get_bucket_df(bucketname, paste0(path, 'BLUP'))
file_list_blup <- dfBucket_blup$Key[grep('heritability', dfBucket_blup$Key)]

file1 <- "ycao1/ITG_data/SweetCorn/BLUP/Fresh_heritability.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/Fresh_heritability_deep.csv"
heritabilityVis(file1, file2, 'Fresh: Heritability Comparison')

file1 <- "ycao1/ITG_data/SweetCorn/BLUP/Proc_heritability.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/Proc_heritability_deep.csv"
heritabilityVis(file1, file2, 'Process: Heritability Comparison')
