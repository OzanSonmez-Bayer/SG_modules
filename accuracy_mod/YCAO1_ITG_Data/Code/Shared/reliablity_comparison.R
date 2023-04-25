# Compare reliabilities between using old and new A matrix

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
pgca_pcm0_plots <- lapply(tomato_trt, FUN = function(x){CompareReliabilityPlot(file1, file2,trait_name = x)})
for (i in 1: length(pgca_pcm0_plots)){
  png(paste0("ProcTomato/plots/Reliability/PCM0/GCA/re_", tomato_trt[i], "_pcm0.png"))
  print(pgca_pcm0_plots[[i]])
  dev.off()
}

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/PGCA_PCM1_2014_2019.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/PGCA_PCM1_2014_2019_deep.csv"
pgca_pcm1_plots <- lapply(tomato_trt, FUN = function(x){CompareReliabilityPlot(file1, file2,trait_name = x)})

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/PSCA_PCM0_2014_2019.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/PSCA_PCM0_2014_2019_deep.csv"
psca_pcm0_plots <- lapply(tomato_trt, FUN = function(x){CompareReliabilityPlot(file1, file2,trait_name = x)})
for (i in 1: length(psca_pcm0_plots)){
  png(paste0("ProcTomato/plots/Reliability/PCM0/SCA/re_", tomato_trt[i], "_pcm0.png"))
  print(psca_pcm0_plots[[i]])
  dev.off()
}

file1 <- "ycao1/ITG_data/ProcTomato/BLUP/PSCA_PCM1_2014_2019.csv"
file2 <- "ycao1/ITG_data/ProcTomato/BLUP/PSCA_PCM1_2014_2019_deep.csv"
psca_pcm1_plots <- lapply(tomato_trt, FUN = function(x){CompareReliabilityPlot(file1, file2,trait_name = x)})


for (i in 1: length(psca_pcm1_plots)){
  png(paste0("ProcTomato/plots/SCA/re_", tomato_trt[i], "_pcm1.png"))
  print(psca_pcm1_plots[[i]])
  dev.off()
}
############################# Sweet Corn ############################
path <- ('ycao1/ITG_data/SweetCorn/')
bucketname <- 'genome-analytics-perm-space'

dfBucket_blup <- get_bucket_df(bucketname, paste0(path, 'BLUP'))
file_list_blup <- dfBucket_blup$Key[c(grep('PGCA', dfBucket_blup$Key), grep('PSCA', dfBucket_blup$Key))]

file_df <- data.frame(file1 = file_list_blup[-grep('_deep', file_list_blup)],
                      file2 = file_list_blup[grep('_deep', file_list_blup)])
file_df

corn_trt <- c('SC_HE','SC_TF','QUAL','SC_PR','SC_SE','HSC','S50D' ,'RTLR' ,'ELEN' ,'EDIA' ,'SDV')

# for (i in 1:nrow(file_df)){
#   pgca_pcm0_plots <- lapply(corn_trt, FUN = function(x){CompareReliabilityPlot(file_df[i, 1], 
#                                                                                file_df[i, 2],
#                                                                                trait_name = x)})
#   for (i in 1: length(pgca_pcm0_plots)){
#     if (grep('GCA_Fresh', file_df[i,1])){
#       png(paste0("Corn/plots/reliability/Fresh/GCA/re_", corn_trt[i], "_fresh.png"))
#     }else if (grep('GCA_Proc', file_df[i,1])){
#       png(paste0("Corn/plots/reliability/Proc/GCA/re_", corn_trt[i], "_proc.png"))
#     }else if (grep('SCA_Fresh', file_df[i,1])){
#       png(paste0("Corn/plots/reliability/Fresh/SCA/re_", corn_trt[i], "_fresh.png"))
#     }else if (grep('SCA_Proc', file_df[i,1])){
#       png(paste0("Corn/plots/reliability/Proc/SCA/re_", corn_trt[i], "_proc.png"))
#     }
#     print(pgca_pcm0_plots[[i]])
#     dev.off()
#   }
#   
#   other_relatives <- lapply(corn_trt, FUN = function(x){CompareReliabilityPlot_other(file_df[i, 1], 
#                                                                                      file_df[i, 2], 
#                                                                                      trait_name = x)})
#   
#   for (i in 1:length(other_relatives)){
#     if (grep('GCA_Fresh', file_df[i,1])){
#       png(paste0("Corn/plots/reliability/Fresh/GCA/re_", corn_trt[i], "_fresh_other.png"))
#     }else if (grep('GCA_Proc', file_df[i,1])){
#       png(paste0("Corn/plots/reliability/Proc/GCA/re_", corn_trt[i], "_proc_other.png"))
#     }else if (grep('SCA_Fresh', file_df[i,1])){
#       png(paste0("Corn/plots/reliability/Fresh/SCA/re_", corn_trt[i], "_fresh_other.png"))
#     }else if (grep('SCA_Proc', file_df[i,1])){
#       png(paste0("Corn/plots/reliability/Proc/SCA/re_", corn_trt[i], "_proc_other.png"))
#     }
#     print(other_relatives[[i]])
#     dev.off()
#   }
#   
# }


file1 <- "ycao1/ITG_data/SweetCorn/BLUP/PGCA_Fresh_2014_2018.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/PGCA_Fresh_2014_2018_deep.csv"
pgca_pcm0_plots <- lapply(corn_trt, FUN = function(x){CompareReliabilityPlot(file1, file2,trait_name = x)})
for (i in 1: length(pgca_pcm0_plots)){
  png(paste0("Corn/plots/reliability/Fresh/GCA/re_", corn_trt[i], "_fresh.png"))
  print(pgca_pcm0_plots[[i]])
  dev.off()
}

other_relatives <- lapply(corn_trt, FUN = function(x){CompareReliabilityPlot_other(file1, file2, trait_name = x)})

for (i in 1:length(other_relatives)){
  png(paste0("Corn/plots/reliability/Fresh/GCA/re_", corn_trt[i], "_fresh_other.png"))
  print(other_relatives[[i]])
  dev.off()
}


file1 <- "ycao1/ITG_data/SweetCorn/BLUP/PGCA_Proc_2014_2018.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/PGCA_Proc_2014_2018_deep.csv"
pgca_pcm1_plots <- lapply(corn_trt, FUN = function(x){CompareReliabilityPlot(file1, file2,trait_name = x)})
for (i in 1: length(pgca_pcm1_plots)){
  png(paste0("Corn/plots/reliability/Proc/GCA/re_", corn_trt[i], "_proc.png"))
  print(pgca_pcm1_plots[[i]])
  dev.off()
}

other_relatives <- lapply(corn_trt, FUN = function(x){CompareReliabilityPlot_other(file1, file2, trait_name = x)})

for (i in 1:length(other_relatives)){
  png(paste0("Corn/plots/reliability/Proc/GCA/re_", corn_trt[i], "_proc_other.png"))
  print(other_relatives[[i]])
  dev.off()
}


file1 <- "ycao1/ITG_data/SweetCorn/BLUP/PSCA_Fresh_2014_2018.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/PSCA_Fresh_2014_2018_deep.csv"
psca_pcm0_plots <- lapply(corn_trt, FUN = function(x){CompareReliabilityPlot(file1, file2,trait_name = x)})
for (i in 1: length(psca_pcm0_plots)){
  png(paste0("Corn/plots/reliability/Fresh/SCA/re_", corn_trt[i], "_fresh.png"))
  print(psca_pcm0_plots[[i]])
  dev.off()
}


file1 <- "ycao1/ITG_data/SweetCorn/BLUP/PSCA_Proc_2014_2018.csv"
file2 <- "ycao1/ITG_data/SweetCorn/BLUP/PSCA_Proc_2014_2018_deep.csv"
psca_pcm1_plots <- lapply(corn_trt, FUN = function(x){CompareReliabilityPlot(file1, file2,trait_name = x)})
for (i in 1: length(psca_pcm1_plots)){
  png(paste0("Corn/plots/reliability/Proc/SCA/re_", corn_trt[i], "_proc.png"))
  print(psca_pcm1_plots[[i]])
  dev.off()
}
