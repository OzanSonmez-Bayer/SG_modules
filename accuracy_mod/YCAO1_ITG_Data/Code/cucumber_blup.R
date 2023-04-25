# Pull clean cucumber data 

# Double check data quality
# raw_data[which(raw_data$PLOT_BID == 'P00000000879547539919167' & raw_data$GROWSEASON == '2015:08' & raw_data$OBSRVTN_REF_CD == 'NETWT'),
#          c('TRAIT_VALUE', 'OBS_DATE_RECORDED', 'OBS_DATE_MODIFIED', 'CROSS_NAME', 'TEST_SET_NAME')]


library(readxl)
library(tidyverse)
library(aws.s3)
library(plyr)
library(dplyr)
library(asreml)
library(asremlPlus)
library(doBy)
library(jsonlite)
library(ggrepel)


if (grepl('WINDOWS', Sys.getenv('windir'))){
  setwd('C:/Users/YCAO1/Desktop/P-1/VEG_ITG/Cucumber')
}

source('Code/Shared/misFunctions.R')
set_df <- read.csv('Cucumber/set_df.csv')
traits <- read_excel('Cucumber/Y3Traits.xlsx')

traits_pcm1 <- traits$FtsTrait[which(traits$PCM1 == 1)]
traits_pcm2 <- traits$FtsTrait[which(traits$PCM2 == 1)]

################################### Pull Data #####################################
crop <- 'Cucumber'
selectedProg <- 'Y3'

raw_data <- data.frame()
years <- unique(set_df$season)[!is.na(unique(set_df$season))]
for (SelectedSeason in years){
  sets <- unique(set_df$set.name[which(set_df$season == SelectedSeason)])
  temp_data <- pullData(season = SelectedSeason, 
                        setlist = sets, 
                        trtlist = traits_pcm1, 
                        crop = crop, 
                        selectedProg = selectedProg)
  
  raw_data <- rbind.fill(raw_data, temp_data)
}

# Check UOM
raw_data <- conv_uom(raw_data)

# Calculated traits
pcm1_dat <- raw_data
pcm1_dat_2 <- data.frame()
for (season in unique(pcm1_dat$GROWSEASON)){
  print(season)
  calc_dat <- pcm1_dat %>% 
    filter(GROWSEASON == season)
  
  calc_dat_1 <-cal_trt(calc_dat)
  pcm1_dat_2 <- rbind(pcm1_dat_2, calc_dat_1)
}

s3save(pcm1_dat_2,bucket = bucket, object = 'ycao1/ITG_data/Cucumber/Raw/pcm1_dat.Rdata') # error
# Pull Data from S3, clean data with calculated traits 
# There are some issues, need check it out
AWS.Authorization('ycao1')
bucket <- 'genome-analytics-perm-space'
infile <- 'ycao1/ITG_data/Cucumber/Raw/pcm1_dat.Rdata'
s3load(object = infile, bucket = bucket)

pcm1_sub <- pcm1_dat %>% 
  filter(DESCRIPTOR %in% c('LATER', 'STEM')) %>% 
  mutate(OBSRVTN_REF_CD = paste0(OBSRVTN_REF_CD,'_',DESCRIPTOR),
         REPETITION = 1)
############################################### Repeated Measures ##############################
Crop <- 'Cucumber'
source('Code/Shared/deepped.R')
source('Code/Shared/BLUPFunctions.R')

unique(pcm1_dat_2$OBSRVTN_REF_CD)

trtList <- c('YNMPA', 'WMFSA', 'FRLGT', 'EXTCO', 'SHAPE', 
             'LNUNI', 'EARLY', 'NODES', 'FRNOD', 'LSHOT', 'AFW_C',
             'NFSA_C', 'WFSA_C')

aggregate_mdl <- c('sum', 'sum', 'tbd', 'tbd', 'tbd', 
                   'tbd', 'single', 'single', 'single', 'single', 'single', 'sum', 'sum')

pcm1_df <- pcm1_dat_2 %>% 
  filter(OBSRVTN_REF_CD %in% trtList)


# Repeated Measures 


# Check repeated measures by grow seasn
# Aggregated (sum or average) them together within each season

pcm1_agg <- data.frame()

for (trt in 1:length(trtList)){
  for (seas in unique(pcm1_df$GROWSEASON)){
    temp_dat <- pcm1_df %>% 
      filter(OBSRVTN_REF_CD == trtList[trt] & GROWSEASON == seas)
    if (length(unique(temp_dat$REPETITION)) > 1){
      temp_agg <- repeatedMeasure(dat = temp_dat, trait =  trtList[trt], season = seas, mdl = aggregate_mdl[trt]) 
      if (nrow(pcm1_agg) == 0){
        pcm1_agg <- temp_agg
      }else{
        pcm1_agg <- rbind.fill(pcm1_agg, temp_agg) 
      }
    }else{
      if (nrow(pcm1_agg) == 0){
        pcm1_agg <- temp_dat
      }else{
        pcm1_agg <- rbind.fill(pcm1_agg, temp_dat) 
      }
    }
  }
}


# analyze aggregated traits together across season, even though it is not the number of repetition

pcm1_agg$OBSRVTN_REF_CD <- gsub('_1.+$', '', pcm1_agg$OBSRVTN_REF_CD)

pcm1_agg_1 <- pcm1_agg
pcm1_agg_1$REPETITION[which(pcm1_agg_1$OBSRVTN_REF_CD %in% c('EARLY', 'NODES', 'FRNOD', 'LSHOT'))] <- 1

############################################### BLUP Calculation ##############################
metric <- 'ABLUP'
AWS.Authorization('ycao1')
bucket <- 'genome-analytics-perm-space'
path <- 'ycao1/ITG_data/Cucumber'
put_folder('/ycao1/ITG_data/Cucumber/BLUP', bucket = bucket)


if (metric == 'ABLUP'){
  pcm1_blup <- Ablup_fun(pcm1_agg_1, n_gen = 5, CropIDs = CropIDs)
  psca_pcm1 <- pcm1_blup[[2]] %>% 
    filter(!is.na(PEDIGREE_NAME))
  pgca_pcm1 <- pcm1_blup[[1]] %>% 
    filter(!is.na(PEDIGREE_NAME))
  
  s3write_using(psca_pcm1, FUN = write.csv,
                bucket = bucket,
                object = 'ycao1/ITG_data/Cucumber/BLUP/PSCA_PCM1_2015_2019_deep.csv')
  
  s3write_using(pgca_pcm1, FUN = write.csv,
                bucket = bucket,
                object = 'ycao1/ITG_data/Cucumber/BLUP/PGCA_PCM1_2015_2019_deep.csv')
  
}else if (metric == 'PBLUP'){
  pcm1_pblup <- pblup_fun(pcm1_agg_1)
  psca_pcm1 <- pcm1_pblup[[2]] %>% 
    filter(!is.na(PEDIGREE_NAME))
  pgca_pcm1 <- pcm1_pblup[[1]] %>% 
    filter(!is.na(PEDIGREE_NAME))
  
  s3write_using(psca_pcm1, FUN = write.csv,
                bucket = bucket,
                object = 'ycao1/ITG_data/Cucumber/BLUP/PSCA_PCM1_2015_2019.csv')
  
  s3write_using(pgca_pcm1, FUN = write.csv,
                bucket = bucket,
                object = 'ycao1/ITG_data/Cucumber/BLUP/PGCA_PCM1_2015_2019.csv')
}



# BLUPs for FRNOD by LATER and STEM

blup_frnod <- pblup_fun(pcm1_sub)
ablup_frnod <- Ablup_fun(pcm1_sub, n_gen = 5, CropIDs = CropIDs)

ablup_frond_stem <- Ablup_fun(pcm1_sub[which(pcm1_sub$OBSRVTN_REF_CD == 'FRNOD_STEM'), ], n_gen = 5, CropIDs = CropIDs)




source('Code/Shared/misFunctions.R')

########################################## SUMMARY STATISTICS ##################################

AWS.Authorization('ycao1')
path <- 'ycao1/ITG_data/Cucumber/'
bucketname <- 'genome-analytics-perm-space'
infile <- paste0('s3://', bucketname, '/', path, 'BLUP/', 'PGCA_PCM1_2015_2019.csv')
# infile <- paste0('s3://', bucketname, '/', path, 'BLUP/', 'h2_pcm1.csv')

csvfile <- get_object(infile)
csvfile_1 <- rawToChar(csvfile)
con <- textConnection(csvfile_1)
pheno <- read.csv(con)
close(con)

psca_pcm1 <- pheno
hr_pcm1_psca <- hrFun(psca_pcm1) 
hr_pcm1_pgca <- hrFun(pgca_pcm1)

s3write_using(hr_pcm1_psca, FUN = write.csv,
              bucket = bucketname,
              object = 'ycao1/ITG_data/Cucumber/BLUP/h2_pcm1.csv')


# Heritability from Deep Pedigree 
hr_pcm1_sca_deep <- hrFun(psca_pcm1)
hr_pcm1_gca_deep <- hrFun(pgca_pcm1)

s3write_using(hr_pcm1_gca_deep, FUN = write.csv,
              bucket = bucketname, 
              object = 'ycao1/ITG_data/Cucumber/BLUP/h2_pcm1_deep.csv')

######################################### COMPARISION ABLUP VS PBLUP ##################################

# AGCA <- A_blup[[1]]
# PGCA <- P_blup[[1]]

AGCA <- pgca_pcm1

PGCA <- pheno[,-1]

ap_merge <- AGCA %>% 
  inner_join(PGCA, by = 'PEDIGREE_NAME')

colnames(ap_merge) <- gsub('.x', '_A',colnames(ap_merge))
colnames(ap_merge) <- gsub('.y', '_P', colnames(ap_merge))


cor(ap_merge$predicted.value_A, ap_merge$predicted.value_P) 


# compare reliability

for (i in 1:length(unique(AGCA$trait))){
  AGCA_sub <- AGCA %>% 
    filter(trait == unique(AGCA$trait)[i])
  
  PGCA_sub <- PGCA %>% 
    filter(trait == unique(PGCA$trait)[i])
  
  ap_reliabity <- rbind(AGCA_sub, PGCA_sub)
  ap_reliabity$mdl <- rep(c('5 Gen','2 Gen'), c(nrow(AGCA_sub), nrow(PGCA_sub)))
  
  ggplot(ap_reliabity, aes(x = mdl, y = reliability, fill = mdl)) + 
    geom_violin() + 
    stat_summary(fun.y = mean, geom = 'point', color = 'blue', size = 5) + 
    xlab('A Matrix') +
    ggtitle(paste(unique(AGCA$trait)[i], ': Reliaibility Comparison'))
  
}


# Heritability Comparison

heritabilityVis(file1 = 'ycao1/ITG_data/Cucumber/BLUP/h2_pcm1.csv',
                file2 = 'ycao1/ITG_data/Cucumber/BLUP/h2_pcm1_deep.csv',
                title_name = 'PCM1: Heritability Comparison')


# Compare BLUPs 
trait_blup <- c('YNMPA_sum', 'WMFSA_sum', 'FRLGT', 'EXTCO', 'SHAPE', 'LNUNI', 'EARLY', 'NODES','FRNOD', 'LSHOT',
                'AFW_C', 'NFSA_C_sum', 'WFSA_C_sum')

file1 <- 'ycao1/ITG_data/Cucumber/BLUP/PSCA_PCM1_2015_2019.csv'
file2 <- 'ycao1/ITG_data/Cucumber/BLUP/PSCA_PCM1_2015_2019_deep.csv'
sca_pcm1_cor <- CompareCorrelation(file1, file2, trait_blup)

file1 <- 'ycao1/ITG_data/Cucumber/BLUP/PGCA_PCM1_2015_2019.csv'
file2 <- 'ycao1/ITG_data/Cucumber/BLUP/PGCA_PCM1_2015_2019_deep.csv'
gca_pcm1_cor <- CompareCorrelation(file1, file2, trait_blup)


pcm1_cor <- rbind(sca_pcm1_cor, gca_pcm1_cor)
pcm1_cor$mdl <- as.factor(rep(c('SCA', 'GCA'), c(nrow(sca_pcm1_cor), nrow(gca_pcm1_cor))))
visualCorr(pcm1_cor, "PCM 1 : BLUP Correlation")



################# Download all the results from AWS #################
AWS.Authorization('ycao1')
path <- 'ycao1/ITG_data/Cucumber/'
bucketname <- 'genome-analytics-perm-space'
infile <- paste0('s3://', bucketname, '/', path, 'BLUP/', 'PGCA_PCM1_2015_2019.csv')
# infile <- paste0('s3://', bucketname, '/', path, 'BLUP/', 'h2_pcm1.csv')

csvfile <- get_object(infile)
csvfile_1 <- rawToChar(csvfile)
con <- textConnection(csvfile_1)
pheno <- read.csv(con)
close(con)


pheno_sub <- pheno %>% 
  filter(trait %in% c('YNMPA_sum', 'WMFSA_sum', 'AFW', 'EARLY', 'NODES', 'LSHOT', 'NFSA_C_sum', 'WFSA_C_sum', 'FRNOD_LATER', 'FRNOD_STEM')) %>% 
  select(-X)

# pheno_sub <- rbind(pheno_sub, ablup_frnod[[2]][which(ablup_frnod[[2]]$trait == 'FRNOD_LATER'), ], ablup_frond_stem[[2]])
# # pheno_sub <- rbind(pheno_sub, blup_frnod[[2]])
# s3write_using(pheno_sub, FUN = write.csv,
#               bucket = bucket,
#               object = 'ycao1/ITG_data/Cucumber/BLUP/PSCA_PCM1_2015_2019.csv')



write.csv(pheno_sub, file = 'PGCA_PCM1_2015_2019.csv', row.names = F)

hr_pcm1_gca_deep <- hrFun(pheno_sub)

write.csv(hr_pcm1_gca_deep, file = 'h2_reliability.csv', row.names = F)


################################### Parallel Computing for Pedigree-based Relationship Matrix #######################

ped_file <- pcm1_df %>% 
  dplyr::select(PEDIGREE_NAME, P1, P2, LTYPE) %>% 
  mutate(pedigree = PEDIGREE_NAME) %>% 
  distinct_all()
dim(ped_file)

new_token(fingerprint = F)
ped_5gen <- make_ped_file(df = ped_file, crop_ids = CropIDs, n_gen = 5)


ped_file_1 <- pcm1_df %>% 
  dplyr::select(PEDIGREE_NAME, P1, P2) %>% 
  rename(Pedigree = PEDIGREE_NAME) %>% 
  mutate(Selfing = 0) %>% 
  distinct_all()
ped_file_1$P1[is.na(ped_file_1$P2)] <- NA
ped_file_1$P2[is.na(ped_file_1$P1)] <- NA


pd.aniv <- asreml.Ainverse(ped_file_1, fgen = c('Selfing',5))
ped_2gen_pblup <- pblup_fun_1(pcm1_agg_1, ainvse_mat = pd.aniv$ginv)

hr_pblup_1 <- hrFun(ped_2gen_pblup[[2]])
hr_pgca_1 <- hrFun(ped_2gen_pblup[[1]])

## Test out repeated measure aggregation 

tst <- pcm1_agg_1 %>% 
  filter(GROWSEASON == '2019:01' & OBSRVTN_REF_CD == 'WFSA_C_sum')

pred_tst <- pblup_fun(tst)


tst_proloog <- pcm1_df[which(pcm1_df$CROSS_NAME == 'PROLOOG' & pcm1_df$GROWSEASON == '2019:01' & pcm1_df$OBSRVTN_REF_CD == 'WFSA_C'),]

tst_proloog_agg <- repeatedMeasure(tst_proloog, 'WFSA_C', '2019:01', 'sum')

# Recalculate 

tst_agg <- pcm1_df %>% 
  filter(OBSRVTN_REF_CD == 'WFSA_C')

tst_agg_1 <- lapply(unique(tst_agg$GROWSEASON), function(x){repeatedMeasure(tst_agg, 'WFSA_C', x, 'sum')})

tst_agg_2 <- do.call(rbind, tst_agg_1)

tst_agg_2_2019 <- tst_agg_2 %>% 
  filter(GROWSEASON == '2019:01')

tst_blup <- pblup_fun(tst_agg_2_2019)


h2h_agg <- read.csv('npheno_agg.csv')
h2h_agg$REPETITION <- 1


tst_blup_1 <- pblup_fun(h2h_agg)
