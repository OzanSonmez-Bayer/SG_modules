# GBLUP for Processing Tomato

library(readxl)
library(tidyverse)
library(ggplot2)
library(aws.s3)
library(jsonlite)
library(Matrix)

library(BGLR)
library(coda)
library(MCMCpack)

source('Code/Shared/BLUPFunctions.R')
source('Code/Shared/genoFunctions.R')
source('Code/Shared/misFunctions.R')
source('Code/Shared/CGS.R')

################# Import Raw Taqman Data #################

procTom_taqman <- read.table('ProcTomato/9Z_TaqManFromCGS_FilteredMarkers.tab', 
                             sep = '\t',
                             fill = TRUE,
                             header = TRUE)
procTom_taqman <- procTom_taqman %>% 
  filter(!is.na(Pedigree))

dim(procTom_taqman)
miss_mrk <- apply(procTom_taqman[,-c(1:2)], MARGIN = 1, FUN = function(x){sum(is.na(x))})
miss_ind <- apply(procTom_taqman[,-c(1:2)], MARGIN = 2, FUN = function(x){sum(is.na(x))})

# remove individuals with 285 missing markers
procTom_taqman_1 <- procTom_taqman[-which(miss_mrk == 287),] # 2607 individuals left
miss_ind <- apply(procTom_taqman_1[,-c(1:2)], MARGIN = 2, FUN = function(x){sum(is.na(x))})
procTom_taqman_2 <- procTom_taqman_1[, -which(miss_ind == 2599)]
geno <- procTom_taqman_2[,-1]


procTom_taqman_2[1:5,1:5]

taqman_full <- procTom_taqman_2 %>% 
  dplyr::select(-Inventory_SampleID) %>% 
  unique()
rownames(taqman_full) <- taqman_full$Pedigree
taqman_full <- taqman_full %>% 
  dplyr::select(-Pedigree)

taqman_full_1 <- gen2Additive(taqman_full)
rownames(taqman_full_1) <- rownames(taqman_full)

G_rrBLUP <- rrBLUP::A.mat((taqman_full_1 - 1), 
                          min.MAF = 0.01,  
                          impute.method = 'mean', 
                          tol = 0.02, 
                          return.imputed = TRUE)

G_imputed <- G_rrBLUP$imputed
hist(diag(G_rrBLUP$A))
range(diag(G_rrBLUP$A)) 



############### Import GCA from Deep Pedigree #############

AWS.Authorization('ycao1')
path <- 'ycao1/ITG_data/ProcTomato/'
bucketname <- 'genome-analytics-perm-space'
infile <- paste0('s3://', bucketname, '/', path, 'BLUP/', 'PGCA_PCM0_2014_2019_deep.csv')
csvfile <- get_object(infile)
csvfile_1 <- rawToChar(csvfile)
con <- textConnection(csvfile_1)
pheno <- read.csv(con)
close(con)


# taqman_full_sub <- taqman_full[which(rownames(taqman_full) %in% as.character(unique(pheno$PEDIGREE_NAME))),]
# taqman_full_sub_012 <- gen2Additive(taqman_full_sub)
# rownames(taqman_full_sub_012) <- rownames(taqman_full_sub)
# 
# 
# G_rrBLUP <- rrBLUP::A.mat((taqman_full_sub_012 - 1), 
#                           min.MAF = 0.01, 
#                           impute.method = 'mean', 
#                           tol = 0.02)
# hist(diag(G_rrBLUP))
# range(diag(G_rrBLUP)) 



############## Match Geno with Pheno #########################

# FBRIX
# MAT, FRFRMH, PLTVG, PLCMP, FQUAL, PYLDPA, SAMPUW, OST, AVJB, LBRIX, FBRIX, AVGHBRIX, FZUNI, MATR, EARLY

trt <- 'FZUNI'
pheno_df <- pheno %>% 
  filter(trait == trt) %>% 
  mutate(dEBV = predicted.value/reliability)

length(intersect(as.character(pheno_df$PEDIGREE_NAME), rownames(G_imputed))) # 496 inbreds have been genotyped

pheno_sub <- pheno_df %>% 
  filter(as.character(PEDIGREE_NAME) %in% rownames(G_imputed)) %>%
  dplyr::select(PEDIGREE_NAME, dEBV)
pheno_sub <- pheno_sub[!duplicated(pheno_sub$PEDIGREE_NAME),]
rownames(pheno_sub) <- pheno_sub$PEDIGREE_NAME
pheno_sub$PEDIGREE_NAME <- as.character(pheno_sub$PEDIGREE_NAME)
# pheno_sub <- pheno_sub %>% 
#   dplyr::select(-PEDIGREE_NAME)


G_sub <- G_imputed[rownames(G_imputed) %in% pheno_sub$PEDIGREE_NAME,]
G_sub <- G_sub[pheno_sub$PEDIGREE_NAME,]


################################## Bayesian LASSO ######################

nIter <- 20000
burnIn <- 2000


# Bayesian LASSO
set.seed(20200422)

ETA <- list(MRK = list(X = G_sub, model = 'BL'))
fmBL <- BGLR(y = pheno_sub$dEBV, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = '')

# plot(abs(fmBL$ETA[[1]]$b),col=4,cex=.5, type='o',main='Bayesian Lasso')
# 
# plot(fmBL$yHat ~ pheno_sub$dEBV, xlab='Observed',ylab='Predicted',col=2,
#      xlim = range(c(fmBL$y, fmBL$yHat)),ylim = range(c(fmBL$y, fmBL$yHat)))
# abline(a=0,b=1,col=4,lwd=2)


r_BL <- cor(pheno_sub$dEBV, fmBL$yHat)
r_BL

# gHat<-G_sub%*%fmBL$ETA[[1]]$b

############################# Bayesian Ridge ###########################

nIter <- 20000
burnIn <- 2000


# Bayesian LASSO
set.seed(20200424)

ETA <- list(MRK = list(X = G_sub, model = 'BRR'))
fmBL <- BGLR(y = pheno_sub$dEBV, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = '')
# 
# plot(abs(fmBL$ETA[[1]]$b),col=4,cex=.5, type='o',main='Bayesian Lasso')
# 
# plot(fmBL$yHat ~ pheno_sub$dEBV, xlab='Observed',ylab='Predicted',col=2,
#      xlim = range(c(fmBL$y, fmBL$yHat)),ylim = range(c(fmBL$y, fmBL$yHat)))
# abline(a=0,b=1,col=4,lwd=2)


r_BL <- cor(pheno_sub$dEBV, fmBL$yHat)
r_BL


############################ Bayes A #################################
nIter <- 20000
burnIn <- 2000


# Bayesian LASSO
set.seed(20200424)

ETA <- list(MRK = list(X = G_sub, model = 'BayesA'))
fmBL <- BGLR(y = pheno_sub$dEBV, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = '')
# 
# plot(abs(fmBL$ETA[[1]]$b),col=4,cex=.5, type='o',main='Bayesian Lasso')
# 
# plot(fmBL$yHat ~ pheno_sub$dEBV, xlab='Observed',ylab='Predicted',col=2,
#      xlim = range(c(fmBL$y, fmBL$yHat)),ylim = range(c(fmBL$y, fmBL$yHat)))
# abline(a=0,b=1,col=4,lwd=2)


r_BL <- cor(pheno_sub$dEBV, fmBL$yHat)
r_BL
