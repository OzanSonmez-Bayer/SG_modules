# Compare BLUP values between 2-gen A matrix and 5-gen A matrix

library(aws.s3)
library(tidyverse)
library(ggplot2)

source('credential.R')
source('misFunctions.R')

############################## Tomato ##################################
path <- ('ycao1/ITG_data/ProcTomato/')
bucketname <- 'genome-analytics-perm-space'

dfBucket_blup <- get_bucket_df(bucketname, paste0(path, 'BLUP'))
file_list_blup <- dfBucket_blup$Key[grep('PCM', dfBucket_blup$Key)]

tomato_trt <- c('MAT','FRFRMH','PLTVG','PLCMP','FQUAL','PYLDPA','SAMPUW',
                'OST','AVJB','LBRIX','FBRIX','AVGHBRIX','FZUNI','MATR','EARLY')

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/PGCA_PCM0_2014_2019.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/PGCA_PCM0_2014_2019_deep.csv"
gca_pcm0_cor <- CompareCorrelation(file1, file2, tomato_trt)
pgca_pcm0_plots <- lapply(tomato_trt, FUN = function(x){CompareBLUP(file1, file2,trait_name = x)})
for (i in 1: length(pgca_pcm0_plots)){
  png(paste0("ProcTomato/plots/BLUP/PCM0/GCA/blup_", tomato_trt[i], "_pcm0.png"))
  print(pgca_pcm0_plots[[i]])
  dev.off()
}

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/PGCA_PCM1_2014_2019.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/PGCA_PCM1_2014_2019_deep.csv"
gca_pcm1_cor <- CompareCorrelation(file1, file2, tomato_trt[which(!tomato_trt %in% c('MATR', "EARLY") )])
pgca_pcm1_plots <- lapply(tomato_trt, FUN = function(x){CompareBLUP(file1, file2,trait_name = x)})
for (i in 1: length(pgca_pcm1_plots)){
  png(paste0("ProcTomato/plots/BLUP/PCM1/GCA/blup_", tomato_trt[i], "_pcm1.png"))
  print(pgca_pcm1_plots[[i]])
  dev.off()
}

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/PSCA_PCM0_2014_2019.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/PSCA_PCM0_2014_2019_deep.csv"
sca_pcm0_cor <- CompareCorrelation(file1, file2, tomato_trt)
psca_pcm0_plots <- lapply(tomato_trt, FUN = function(x){CompareBLUP(file1, file2,trait_name = x)})
for (i in 1: length(psca_pcm0_plots)){
  png(paste0("ProcTomato/plots/BLUP/PCM0/SCA/blup_", tomato_trt[i], "_pcm0.png"))
  print(psca_pcm0_plots[[i]])
  dev.off()
}

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/PSCA_PCM1_2014_2019.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/PSCA_PCM1_2014_2019_deep.csv"
sca_pcm1_cor <- CompareCorrelation(file1, file2, tomato_trt[which(!tomato_trt %in% c('MATR', "EARLY") )])
psca_pcm1_plots <- lapply(tomato_trt, FUN = function(x){CompareBLUP(file1, file2,trait_name = x)})
for (i in 1: length(psca_pcm1_plots)){
  png(paste0("ProcTomato/plots/BLUP/PCM0/SCA/blup_", tomato_trt[i], "_pcm1.png"))
  print(psca_pcm1_plots[[i]])
  dev.off()
}

# Correlation visualization
pcm0_cor <- rbind(sca_pcm0_cor, gca_pcm0_cor)
pcm0_cor$mdl <- as.factor(rep(c('SCA', 'GCA'), c(nrow(sca_pcm0_cor), nrow(gca_pcm0_cor))))
visualCorr(pcm0_cor, "PCM 0 : BLUP Correlation")

pcm1_cor <- rbind(sca_pcm1_cor, gca_pcm1_cor)
pcm1_cor$mdl <- as.factor(rep(c('SCA', 'GCA'), c(nrow(sca_pcm1_cor), nrow(gca_pcm1_cor))))
visualCorr(pcm1_cor)


# Check individuals without direct observational data from hybrids
linesCheck(file1 = file2)

############################# Sweet Corn ############################
path <- ('ycao1/ITG_data/SweetCorn/')
bucketname <- 'genome-analytics-perm-space'

dfBucket_blup <- get_bucket_df(bucketname, paste0(path, 'BLUP'))
file_list_blup <- dfBucket_blup$Key[c(grep('PGCA', dfBucket_blup$Key), grep('PSCA', dfBucket_blup$Key))]

corn_trt <- c('SC_HE','SC_TF','QUAL','SC_PR','SC_SE','HSC','S50D' ,'RTLR' ,'ELEN' ,'EDIA' ,'SDV')

file1 <- "ycao1/ITG_data/SweetCorn/BLUP/PGCA_Fresh_2014_2018.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/PGCA_Fresh_2014_2018_deep.csv"
gca_pcm0_cor <- CompareCorrelation(file1, file2, corn_trt)
pgca_pcm0_plots <- lapply(corn_trt, FUN = function(x){CompareBLUP(file1, file2,trait_name = x)})
for (i in 1: length(pgca_pcm0_plots)){
  png(paste0("Corn/plots/BLUP/Fresh/GCA/blup_", corn_trt[i], "_fresh.png"))
  print(pgca_pcm0_plots[[i]])
  dev.off()
}


file1 <- "ycao1/ITG_data/SweetCorn/BLUP/PGCA_Proc_2014_2018.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/PGCA_Proc_2014_2018_deep.csv"
gca_pcm1_cor <- CompareCorrelation(file1, file2, corn_trt)
pgca_pcm1_plots <- lapply(corn_trt, FUN = function(x){CompareBLUP(file1, file2,trait_name = x)})
for (i in 1: length(pgca_pcm1_plots)){
  png(paste0("Corn/plots/BLUP/Proc/GCA/blup_", corn_trt[i], "_fresh.png"))
  print(pgca_pcm1_plots[[i]])
  dev.off()
}


file1 <- "ycao1/ITG_data/SweetCorn/BLUP/PSCA_Fresh_2014_2018.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/PSCA_Fresh_2014_2018_deep.csv"
sca_pcm0_cor <- CompareCorrelation(file1, file2, corn_trt)
psca_pcm0_plots <- lapply(corn_trt, FUN = function(x){CompareBLUP(file1, file2,trait_name = x)})
for (i in 1: length(psca_pcm0_plots)){
  png(paste0("Corn/plots/BLUP/Fresh/SCA/blup_", corn_trt[i], "_fresh.png"))
  print(psca_pcm0_plots[[i]])
  dev.off()
}

file1 <- "ycao1/ITG_data/SweetCorn/BLUP/PSCA_Proc_2014_2018.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/PSCA_Proc_2014_2018_deep.csv"
sca_pcm1_cor <- CompareCorrelation(file1, file2, corn_trt)
psca_pcm1_plots <- lapply(corn_trt, FUN = function(x){CompareBLUP(file1, file2,trait_name = x)})
for (i in 1: length(psca_pcm1_plots)){
  png(paste0("Corn/plots/BLUP/Proc/SCA/blup_", corn_trt[i], "_proc.png"))
  print(psca_pcm1_plots[[i]])
  dev.off()
}

# Correlation visualization
pcm0_cor <- rbind(sca_pcm0_cor, gca_pcm0_cor)
pcm0_cor$mdl <- as.factor(rep(c('SCA', 'GCA'), c(nrow(sca_pcm0_cor), nrow(gca_pcm0_cor))))
visualCorr(pcm0_cor, "Fresh : BLUP Correlation")

pcm1_cor <- rbind(sca_pcm1_cor, gca_pcm1_cor)
pcm1_cor$mdl <- as.factor(rep(c('SCA', 'GCA'), c(nrow(sca_pcm1_cor), nrow(gca_pcm1_cor))))
visualCorr(pcm1_cor, "Process : BLUP Correlation")
