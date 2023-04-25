# Genetic gain
# Yearly BLUP from 2014:04 to 2019:04
# PYLDPA

library(aws.s3)
library(httr)
library(asreml)
library(asremlPlus)
library(tidyverse)

install.packages('/mnt/packages/doBy_4.6-2.tar.gz', type = 'source', repos = NULL)
library(doBy)


source('Credential/vaultCredentials.R')
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")


pheno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/pheno_training.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


seasons_sets <- list(
  '2014' = list('year' =  c('2014:04'), 'set' = c('AO_PCM2', 'CK_PCM2')), 
  '2015' = list('year' = c('2015:04'), 'set' = c('CKRT', 'AORT')),
  '2016' = list('year' = c('2016:04'), 'set' = c('PCM1')),
  '2017' = list('year' = c('2017:04'), 'set' = c('PCM1')),
  '2018' = list('year' = c('2018:04'), 'set' = c('PCM1')),
  '2019' = list('year' = c('2019:04'), 'set' = c('PCM1', 'PCM1LATE', 'PMC1TI'))
)


# Run BLUP
source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')

avg_gca <- data.frame(year = names(seasons_sets), average_gca = NA, average_pred = NA)

for (i in 1:length(seasons_sets)){
  temp <- pheno %>% 
    # filter(OBSRVTN_REF_CD == 'AVJB') %>%
    filter(OBSRVTN_REF_CD == 'PYLDPA') %>%
    filter(GROWSEASON == seasons_sets[[i]][['year']] & TEST_SET_NAME %in% seasons_sets[[i]][['set']])
  
  blup_temp <- pblup_fun(dataset = temp, crop = 'Tomato')
  
  gca_df <- blup_temp[[1]]
  
  write.csv(gca_df, file = paste0('/mnt/ProcTomato/', names(seasons_sets)[i], '_GCA_PYLDPA.csv'), row.names = F)
  
  avg_gca[i, 'average_gca'] <- mean(gca_df$BLUP, na.rm = T)
  avg_gca[i, 'average_pred'] <- mean(gca_df$predicted.value, na.rm = T)
}

avg_gca
write.csv(avg_gca, file = 'genetic_gain_by_avg_gca_PYLDPA.csv', row.names = F)


library(ggplot2)

ggplot(avg_gca, aes(x = as.numeric(as.character(year)), y = average_gca)) +
  geom_line( color="#69b3a2", size=1, alpha=0.9, linetype=2) +
  geom_point(size = 3) + 
  ggtitle("Average GCA from 2014 to 2019") + 
  theme_classic() + 
  ylab('Average GCA') + 
  xlab('Year')

ggplot(avg_gca, aes(x = as.numeric(as.character(year)), y = average_pred)) +
  geom_line( color="#69b3a2", size=1, alpha=0.9, linetype=2) +
  geom_point(size = 2) + 
  ggtitle("Average Predicted Value from 2014 to 2019") + 
  theme_classic() + 
  ylab('Average Predicted Value (tons/acre)') + 
  xlab('Year')


dat <- read.csv('/mnt/ProcTomato/2014_GCA_AVJB.csv')
length(unique(dat$PEDIGREE_NAME))
hist(dat$predicted.value)

for (i in 1:length(seasons_sets)){
  temp <- pheno %>% 
    filter(OBSRVTN_REF_CD == 'PYLDPA') %>% 
    filter(GROWSEASON == seasons_sets[[i]][['year']] & TEST_SET_NAME %in% seasons_sets[[i]][['set']])
  
 print(seasons_sets[[i]][['year']])
 print(unique(temp$PMC))
}




############################ G X E Model ############################

library(aws.s3)
library(tidyverse)
library(dplyr)
library(httr)
library(plyr)
library(data.table)


bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
crop <- 'Tomato'

# Credentials

source('Credential/vaultCredentials.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

pheno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/pheno_training.rds', 
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

pheno_pcm1 <- pheno %>% 
  filter(TEST_SET_NAME %in% c('PCM1', 'PCM1TI', 'PMC1LATE', 'CKRT', 'AORT', 'AO_PCM2', 'CK_PCM2')) %>% 
  filter(OBSRVTN_REF_CD %in% c('PYLDPA', 'AVJB'))

unique(pheno_pcm1$OBSRVTN_REF_CD)

rm(pheno)
gc()


library(asreml)
library(asremlPlus)

source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')


pheno_pcm1$TREP <- paste(pheno_pcm1$REP_NUMBER,pheno_pcm1$TRACK_ID,sep='_')

bdat <- pheno_pcm1 %>% 
  filter(OBSRVTN_REF_CD == 'AVJB') %>% 
  mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
  mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
  mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME)) %>% 
  mutate(BR_FIELD_ID = as.character(BR_FIELD_ID)) %>% 
  mutate(ORIGIN = as.character(ORIGIN))


if(sum(bdat$ORIGIN == '')>0){
  bdat$ORIGIN[which(bdat$ORIGIN == '')] <- bdat$PEDIGREE_NAME[which(bdat$ORIGIN == '')]
}

if (!'P1' %in% colnames(bdat)){
  print("parentage pedigree doesn't exit")
  bdat <- bdat %>% 
    separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
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
  
}else{
  print("parentage pedigree exist")
  colnames(bdat)[which(colnames(bdat) == 'parent.female_pedigree_name')] = 'P1'
  colnames(bdat)[which(colnames(bdat) == 'parent.male_pedigree_name')] = 'P2'
  
}
if (sum(is.na(bdat$P1)) >0) {
  bdat[is.na(bdat$P1), 'P2'] <- NA
}

if (sum(is.na(bdat$P2)) >0 ){
  bdat[is.na(bdat$P2), 'P1'] <- NA
}

bdat$P1 <- as.character(bdat$P1)
bdat$P2 <- as.character(bdat$P2)


pd <- unique(data.frame(bdat$PEDIGREE_NAME, bdat$P1, bdat$P2))
pd <- pd[!duplicated(pd[,1]),]
names(pd) <- c('Pedigree', 'P1', 'P2')
pd$Selfing <- rep(0, rep = nrow(pd))

pd.aniv=asreml.Ainverse(pd,fgen = c('Selfing',5))


bdat <- bdat %>% 
  mutate(PEDIGREE_NAME = as.factor(PEDIGREE_NAME))


mdl_1 <- asreml(TRAIT_VALUE ~ 1, 
                random = ~ fa(GROWSEASON,2):ped(PEDIGREE_NAME) + BR_FIELD_ID:TREP, 
                na.method.X = 'include',
                trace = FALSE,
                ginverse = list(PEDIGREE_NAME = pd.aniv$ginv),
                workspace = 320e+6, 
                pworkspace = 320e+6,
                data = bdat)

mdl_1 <- update(mdl_1)

model_id <- mdl_1

grand_mean <- data.frame(coef(model_id)$fixed)
pred_ebv_1 <-data.frame(summary(model_id,all=T)$coef.rand)
pred_ebv_1$predicted.value <- as.numeric(round(pred_ebv_1$solution,4))
pred_ebv_1$predicted.value <- pred_ebv_1$predicted.value + grand_mean$effect
pred_ebv_1$standard.error <- as.numeric(round(pred_ebv_1$std.error,4))

head(pred_ebv_1)

prediction_mdl1 <- pred_ebv_1[grep('fa|at', rownames(pred_ebv_1)),]
prediction_mdl1$CROSS_NAME <- sub(".*NAME)_", "", rownames(prediction_mdl1))
prediction_mdl1$GROWSEASON <- unlist(lapply(rownames(prediction_mdl1), 
                                            function(x){unlist(strsplit(x, split = '_'))[2]}))
prediction_mdl1$GROWSEASON <- unlist(lapply(prediction_mdl1$GROWSEASON, 
                                            function(x){unlist(strsplit(x, split = ':ped'))[1]}))
# prediction_mdl1$Location <- unlist(lapply(prediction_mdl1$Location, 
#                                           function(x){unlist(strsplit(x, split = '_'))[1]}))
# prediction_mdl1$Location <- gsub('field_chk_|_N', '', prediction_mdl1$Location)
# 

prediction_mdl1 <- prediction_mdl1 %>% 
  select(CROSS_NAME, GROWSEASON, predicted.value, solution, standard.error)
rownames(prediction_mdl1) <- NULL

head(prediction_mdl1)

prediction_mdl1$mdl <- NA
prediction_mdl1$mdl[prediction_mdl1$CROSS_NAME %in% as.character(unique(bdat$PEDIGREE_NAME))] <- 'SCA'
prediction_mdl1$mdl[!prediction_mdl1$CROSS_NAME %in% as.character(unique(bdat$PEDIGREE_NAME))] <- 'GCA'

head(prediction_mdl1)

# Average breeding value by GROWSEASON

overall_mean <- mean(bdat$TRAIT_VALUE, na.rm = T)

PYLDPA_gg <- prediction_mdl1 %>% 
  filter(mdl == 'GCA') %>% 
  mutate(GROWSEASON = as.factor(GROWSEASON)) %>% 
  group_by(GROWSEASON) %>% 
  dplyr::summarise(avg_ebv_PYLDPA = mean(solution, na.rm = T)) %>% 
  mutate(avg_ebv_PYLDPA_1 = avg_ebv_PYLDPA + overall_mean) %>% 
  ungroup() %>% 
  filter(!GROWSEASON %in% c('Comp1', 'Comp2')) %>% 
  as.data.frame()

head(PYLDPA_gg)

PYLDPA_gg <- PYLDPA_gg %>% 
  mutate(YEAR = as.numeric(str_sub(GROWSEASON, start = 1L, end = 4L)))

library(ggplot2)

ggplot(PYLDPA_gg, aes(x = YEAR, y = avg_ebv_PYLDPA_1)) +
  geom_line(color="#69b3a2", size=1, alpha=0.9, linetype=2) +
  geom_point(size = 2) + 
  ggtitle("Average Predicted Value of PYLDPA from 2014 to 2019") + 
  theme_classic() + 
  ylab('PYLDPA (tons/acre)') + 
  xlab('Year')

write.csv(prediction_mdl1, file = '/mnt/ProcTomato/PYLDPA_GCA_GXEModel_06062021.csv', row.names = F)
write.csv(PYLDPA_gg, file = '/mnt/ProcTomato/PYLDPA_GG_GXEModel_06062021.csv', row.names = F)

######## AVJB ##########


AVJB_gg <- prediction_mdl1 %>% 
  filter(mdl == 'GCA') %>% 
  mutate(GROWSEASON = as.factor(GROWSEASON)) %>% 
  group_by(GROWSEASON) %>% 
  dplyr::summarise(avg_ebv_AVJB = mean(solution, na.rm = T)) %>% 
  mutate(avg_ebv_AVJB_1 = avg_ebv_AVJB + overall_mean) %>% 
  ungroup() %>% 
  filter(!GROWSEASON %in% c('Comp1', 'Comp2')) %>% 
  as.data.frame()

head(AVJB_gg)

AVJB_gg <- AVJB_gg %>% 
  mutate(YEAR = as.numeric(str_sub(GROWSEASON, start = 1L, end = 4L)))

library(ggplot2)

ggplot(AVJB_gg, aes(x = YEAR, y = avg_ebv_AVJB_1)) +
  geom_line(color="#69b3a2", size=1, alpha=0.9, linetype=2) +
  geom_point(size = 2) + 
  ggtitle("Average Predicted Value of AVJB from 2014 to 2019") + 
  theme_classic() + 
  ylab('AVJB') + 
  xlab('Year')

write.csv(prediction_mdl1, file = '/mnt/ProcTomato/AVJB_GCA_GXEModel_06062021.csv', row.names = F)
write.csv(AVJB_gg,file = '/mnt/ProcTomato/AVJB_GG_GXEModel_06062021.csv', row.names = F)

