# GBLUP

# Tomato

library(rrBLUP)
library(tidyverse)
library(plyr)
library(aws.s3)


source('genoFunctions.R')

#################################### Phenotypic Data ##############################

source('credential.R')

bucektname <- 'genome-analytics-perm-space'
path <- 'ycao1/ITG_data/ProcTomato/BLUP/'

files <- get_bucket(bucket = bucektname, prefix = path)

# pheno data is GCA_PCM0_2014_2018
pheno_dat <- read.csv(text = rawToChar(get_object(files[[2]]$Key, bucket = bucektname)))

deregress_dat <- pheno_dat %>% 
  select(PEDIGREE_NAME, BLUP, reliability, h2, trait) %>% 
  mutate(dEBV = BLUP/reliability)

write.csv(deregress_dat, 'dEBV.csv')
################################### Genotypic Data #################################

tst <- read.csv('dEBV.csv')
paste0(unique(tst$PEDIGREE_NAME)[801:1200], collapse = ',')
