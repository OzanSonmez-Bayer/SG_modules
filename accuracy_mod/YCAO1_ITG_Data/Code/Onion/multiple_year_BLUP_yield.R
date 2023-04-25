# This script is to generate multiple year report for mid-parent values
# It should go through AOB 

library(data.table)
library(reshape2)
library(dplyr)
library(DT)
library(sparkline)
library(car)
library(outliers)
library(GGally)
library(tidyr)
library(asreml)
library(asremlPlus)
library(measurements)
library(abind)
library(arm)
library(aws.s3)
library(stringr)
library(RColorBrewer)
library(R.utils)
library(parallel)
library(httr)

install.packages('/mnt/packages/doBy_4.6-2.tar.gz', repos = NULL, type = 'source')

library(doBy)

source('/mnt/Credential/vaultCredentials.R')
source('/mnt/Credential/Credential_sdi.R')
source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')
# source('/repos/ITG_Data/Code/Shared/deepped.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

Sys.setenv("CLIENT_ID" = CLIENT_ID,
           "CLIENT_SECRET" = CLIENT_SECRET)


bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
Crop <- 'Onion'
program <- c('UH', 'RV')

source('/repos/ITG_Data/Code/Shared/deepped.R')

years <- c('2017','2018','2019', '2020', '2021')
AllData <- data.frame()
for (selectedYear in years){
  infile <- paste0(path,Crop,"/",selectedYear,"/AllYearData.RData")
  s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")
  AllData <- rbind(AllData,YearPhenos)
}

# seasons <- c('2016:04', '2016:07',
#              '2017:04', '2017:07', 
#              '2018:04', '2018:07', 
#              '2019:04', '2019:07', 
#              '2020:04', '2020:07')

trial_type <- 'YieldTrial'

AllData_1 <- AllData %>% 
  filter(SET_OWNER_PROG %in% program & TRIAL_TYPE == trial_type)

PMC <- list(
  'Yellow' = c('OMY-M','OMY', 'OSG', 'OSY-M-EARLY', 'OSY', 'OSY-M-MAIN'),
  'Red' = c('OMR-M', 'OMR', 'OSR', 'OSR-M'),
  'White' = c('OMW-M', 'OMW', 'OSW', 'OSW-M-EARLY')
)

PMC <- c(PMC[['Yellow']], PMC[['Red']], PMC[['White']])
trait <- c('BLSNC', 'SCTOT', 'TBLNO', 'TBLPCT', 'BSHPE')

# trait <- list(
#   'Yellow' = c('BEXYC', 'BFIRM', 'BFSCQ', 'BLSNC', 'BSHPE', 'BSUNI', 'LEREC', 'MUNI', 'PLTVG', 'PYRUV', 'QBRIX', 'TBLNO',
#                'TBLPCT', 'THT', 'TLCLR', 'BUMAT'),
#   'Red' = c('BEXRC', 'BFIRM', 'BFSCQ', 'BINRC', 'BLSNC', 'BSHPE', 'BSUNI', 'LEREC', 'MUNI', 'PLTVG', 'PYRUV', 'QBRIX', 'TBLNO',
#             'TBLPCT', 'THT', 'TLCLR', 'BUMAT'), 
#   'White' = c('BFIRM', 'BFSCQ', 'BLSNC', 'BSHPE', 'BSUNI', 'LEREC', 'MUNI', 'PLTVG', 'PYRUV', 'QBRIX', 'TBLNO',
#               'TBLPCT', 'THT', 'TLCLR', 'BWHCL', 'BWHGR', 'BUMAT')
# )

# market <- 'Yellow'

dat <- AllData_1 %>% 
  filter(OBSRVTN_REF_CD %in% trait) %>% 
  filter(SUB_COUNTRY %in% c('Texas', 'California'))
dat <- dat[grep(paste(PMC, collapse = '|'), dat$PMC),]

# Create a new trait, %single center (BLSNC/SCTOT)

blsnc <- dat %>% 
  filter(OBSRVTN_REF_CD == 'BLSNC') %>% 
  rename(BLSNC = TRAIT_VALUE)

sctot <- dat %>% 
  filter(OBSRVTN_REF_CD == 'SCTOT') %>% 
  rename(SCTOT = TRAIT_VALUE)

blsnc_sctot <- blsnc %>% 
  inner_join(sctot, by = c('PLOT_ID', 'GERMPLASM_ID')) %>% 
  mutate(TRAIT_VALUE = as.numeric(BLSNC)/as.numeric(SCTOT) * 100) %>% 
  mutate(BLSNC = TRAIT_VALUE)


blsnc_sctot <- blsnc_sctot[,c(1:ncol(blsnc))]
  

blsnc_sctot$OBSRVTN_REF_CD.x <- 'BLSNC/SCTOT'


colnames(blsnc_sctot) <- colnames(blsnc)
colnames(blsnc_sctot)[colnames(blsnc_sctot) == 'BLSNC'] <- 'TRAIT_VALUE'

dat <- rbind(dat, blsnc_sctot)

# tst <- Ablup_fun(dat[dat$OBSRVTN_REF_CD == 'BEXYC',], 5, CropIDs, crop = Crop)

BLUP <- Ablup_fun(dat, 5, CropIDs, crop = Crop)
GCA_df <- BLUP[[1]]
SCA_df <- BLUP[[2]]

# The same pedgiree may get multiple BLUPs. Take average for now? Why?

# Add origin 
s3load(object = 'Onion/Onion_IDs.RData', bucket = 'veg-apd-sdi-predictiveanalytcs-prod-reference-data')

GCA_df_origin <- GCA_df %>% 
  left_join(CropIDs[, c('M.GERMPLASM.PEDIGREE', 'M.GERMPLASM.ORIGIN')], by = c('PEDIGREE_NAME' = 'M.GERMPLASM.PEDIGREE')) %>% 
  unique()

SCA_df_origin <- SCA_df %>% 
  left_join(CropIDs[, c('M.GERMPLASM.PEDIGREE', 'M.GERMPLASM.ORIGIN')], by = c('PEDIGREE_NAME' = 'M.GERMPLASM.PEDIGREE')) %>% 
  unique()

s3write_using(FUN = write.csv, x = GCA_df_origin, 
              object = 'ycao1/OriginSelection/Onion/MultipleYear_BLUP/AllMCT_GCA_longformat.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

s3write_using(FUN = write.csv, x = SCA_df_origin, 
              object = 'ycao1/OriginSelection/Onion/MultipleYear_BLUP/AllMCT_SCA_longformat.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

# Wide format 
# 
# GCA_df <- s3read_using(FUN = read.csv,
#                        object = 'ycao1/OriginSelection/Onion/MultipleYear_BLUP/Red_GCA_longformat.csv',
#                        bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
# 
# 
# GCA_df <- GCA_df[,-1]
# 
# head(GCA_df)

# GCA_df <- GCA_df %>% 
#   dplyr::select(PEDIGREE_NAME, predicted.value, BLUP, reliability,trait, n, h2) %>% 
#   mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME))

# The same pedigree has different GCA values, why? Take average for now. 

GCA_df <- GCA_df %>% 
  dplyr::select(PEDIGREE_NAME, predicted.value, trait) %>% 
  mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME)) %>% 
  group_by(trait, PEDIGREE_NAME) %>% 
  dplyr::summarise(predicted.value = mean(predicted.value)) %>% 
  ungroup()

gca_df_wide <- data.frame()

for (tt in unique(dat$OBSRVTN_REF_CD)){
  if (tt %in% GCA_df$trait){
    print(tt)
    temp <- GCA_df %>% 
      filter(trait == tt) %>% 
      dplyr::select(-trait)
    
    colnames(temp)[which(colnames(temp) != 'PEDIGREE_NAME')] <- paste(tt, colnames(temp)[which(colnames(temp) != 'PEDIGREE_NAME')], sep = '_')
    
    if (nrow(gca_df_wide) == 0){
      gca_df_wide <- temp
    }else{
      gca_df_wide <- gca_df_wide %>% 
        full_join(temp, by = 'PEDIGREE_NAME') %>% 
        unique()
    }
    
    rm(temp)
    gc()
  }
}

gca_df_wide <- gca_df_wide %>% 
  left_join(CropIDs[, c('M.GERMPLASM.PEDIGREE', 'M.GERMPLASM.ORIGIN')], by = c('PEDIGREE_NAME' = 'M.GERMPLASM.PEDIGREE')) %>% 
  unique()

gca_df_wide <- gca_df_wide[!duplicated(gca_df_wide$PEDIGREE_NAME),]


s3write_using(FUN = write.csv, x = gca_df_wide, 
              object = 'ycao1/OriginSelection/Onion/MultipleYear_BLUP/AllMCT_GCA_wideformat.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


################################ Mid-parent Value ########################

# Generate Mid-parent value 
install.packages('arrangements')
library(arrangements)

midparent_fun <- function(gca_df_wide, market, CropIDs){
  
  mid_parent <- data.frame()
  
  for (trt in trait[[market]]){
    
    if (sum(grepl(trt, colnames(gca_df_wide))) == 1){
      
      print(trt)
      
      xcol <- paste0(trt, '_predicted.value')
      # temp <- gca_df_wide %>% 
      #   dplyr::select(PEDIGREE_NAME, xcol)
      # 
      temp <- gca_df_wide[, c('PEDIGREE_NAME', xcol)]
      temp <- unique(temp)
      
      comb_ped <- combinations(unique(as.character(temp$PEDIGREE_NAME)),2)
      dim(comb_ped)
      
      colnames(comb_ped) <- c('P1_PEDIGREE', 'P2_PEDIGREE')
      
      comb_ped <- as.data.frame(comb_ped)
      comb_ped$P1_PEDIGREE <- as.character(comb_ped$P1_PEDIGREE)
      comb_ped$P2_PEDIGREE <- as.character(comb_ped$P2_PEDIGREE)
      
      temp_1 <- temp %>% 
        right_join(comb_ped, by = c('PEDIGREE_NAME' = 'P1_PEDIGREE')) %>% 
        unique()
      
      colnames(temp_1)[which(colnames(temp_1) == xcol)] <- paste0('P1_', trt)
      
      temp_2 <- temp_1 %>% 
        right_join(comb_ped, by = 'P2_PEDIGREE') %>% 
        dplyr::select(-P1_PEDIGREE) %>% 
        left_join(temp, by = c('P2_PEDIGREE' = 'PEDIGREE_NAME')) %>% 
        unique() %>% 
        rename(P1_PEDIGREE = PEDIGREE_NAME)
      
      colnames(temp_2)[which(colnames(temp_2) == xcol)] <- paste0('P2_', trt)
      
      temp_2[paste0(trt,'_midparent')] <- (temp_2[,paste0('P1_', trt)] + temp_2[, paste0('P2_', trt)])/2
      
      print(dim(temp_2))
      
      if (nrow(mid_parent) == 0){
        mid_parent = temp_2
      }else{
        mid_parent = mid_parent %>% 
          full_join(temp_2, by = c('P1_PEDIGREE', 'P2_PEDIGREE'))
      }
      rm(temp)
      rm(temp_1)
      rm(temp_2)
      gc() 
    }
  }
  
  mid_parent_1 <- mid_parent %>% 
    mutate(clean_origin = paste0(P1_PEDIGREE, '__', P2_PEDIGREE)) %>% 
    left_join(CropIDs[, c('M.GERMPLASM.PEDIGREE', 'M.GERMPLASM.ORIGIN')], by = c('P1_PEDIGREE' = 'M.GERMPLASM.PEDIGREE')) %>% 
    dplyr::rename(P1_Origin = M.GERMPLASM.ORIGIN) %>% 
    unique() %>% 
    left_join(CropIDs[, c('M.GERMPLASM.PEDIGREE', 'M.GERMPLASM.ORIGIN')], by = c('P2_PEDIGREE' = 'M.GERMPLASM.PEDIGREE')) %>% 
    dplyr::rename(P2_Origin = M.GERMPLASM.ORIGIN) %>% 
    unique() %>% 
    dplyr::select(clean_origin, P1_PEDIGREE, P2_PEDIGREE, P1_Origin, P2_Origin, 
                  colnames(mid_parent)[which(!colnames(mid_parent) %in% c('P1_PEDIGREE', 'P2_PEDIGREE'))]) 
  return(mid_parent_1)
}


check_existing <- function(sca_df, midparent_df){
  sca_df$clean_origin = paste
  
  midparent_df_1 <- midparent_df %>% 
    mutate(existing_status = NA) %>% 
    mutate(existing_status = if_else(as.character(P1_PEDIGREE) %in% as.character(sca_df$PEDIGREE_NAME), 1, 0))
  # mutate(existing_status = if_else(as.character(clean_origin) %in% as.character(sca_df$PEDIGREE_NAME), 1, 0))
  
  return(midparent_df_1)
}


trait <- list(
  'Yellow' = c('BEXYC', 'BFIRM', 'BFSCQ', 'BLSNC', 'BSHPE', 'BSUNI', 'LEREC', 'MUNI', 'PLTVG', 'PYRUV', 'QBRIX', 'TBLNO',
               'TBLPCT', 'THT', 'TLCLR', 'BUMAT'),
  'Red' = c('BEXRC', 'BFIRM', 'BFSCQ', 'BINRC', 'BLSNC', 'BSHPE', 'BSUNI', 'LEREC', 'MUNI', 'PLTVG', 'PYRUV', 'QBRIX', 'TBLNO',
            'TBLPCT', 'THT', 'TLCLR', 'BUMAT'), 
  'White' = c('BFIRM', 'BFSCQ', 'BLSNC', 'BSHPE', 'BSUNI', 'LEREC', 'MUNI', 'PLTVG', 'PYRUV', 'QBRIX', 'TBLNO',
              'TBLPCT', 'THT', 'TLCLR', 'BWHCL', 'BWHGR', 'BUMAT')
)



gca_df_wide <- s3read_using(FUN = read.csv, 
                            object = 'ycao1/OriginSelection/Onion/MultipleYear_BLUP/White_GCA_wideformat.csv', 
                            bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

sca_df <- s3read_using(FUN = read.csv, 
                       object = 'ycao1/OriginSelection/Onion/MultipleYear_BLUP/White_SCA_longformat.csv', 
                       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

add_p1p2 <- function(bdat){
  bdat <- bdat %>% 
    separate(M.GERMPLASM.ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
  p1 <- unlist(lapply(bdat$P1, FUN = pedSingle))
  p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  bdat$P1 <- p1
  bdat$P2[bdat$P2=='']=NA
  p2 <- unlist(lapply(bdat$P2, FUN = pedSingle))
  p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
  bdat$P2 <- p2
  
  return(bdat)
}

sca_df <- add_p1p2(sca_df)

market <- 'White'

devX_table <- midparent_fun(gca_df_wide = gca_df_wide, market = market, CropIDs = CropIDs)

mid_parent_1 <- check_existing(sca_df = sca_df, midparent_df = devX_table)


s3write_using(FUN = write.csv, x = mid_parent_1, 
              object = 'ycao1/OriginSelection/Onion/DiallelTable/White_devX.csv', 
              bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')



