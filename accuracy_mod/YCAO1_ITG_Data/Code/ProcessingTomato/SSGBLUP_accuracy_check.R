# Check accuracy 
library(dplyr)
library(aws.s3)
library(httr)
library(ggplot2)
library(readxl)

source('/repos/BLUPF90/BLUPF90/funcStore_deep_ITG.R')

loadcrential()
Alldata <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/pheno_training.rds',
                     bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

geno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/geno_all.rds',
                  bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


trait <- 'SAMPUW'
ssgblup_path <- '/mnt/Tomato/SSGBLUP'
ablup_path <- '/mnt/Tomato/ABLUP'


observed_pheno <- Alldata %>% 
  filter(GROWSEASON == '2019:04' & OBSRVTN_REF_CD == trait) %>% 
  select(PEDIGREE_NAME, TRAIT_VALUE) %>% 
  mutate(TRAIT_VALUE  = as.numeric(TRAIT_VALUE))

load(paste0(ssgblup_path, '/', trait, '_deep_fold_1', '/', 'BLUPF90_removeYearNone.RData'))

load(paste0(ablup_path, '/', trait, '_deep_fold_1', '/', 'BLUPF90_removeYearNone.RData'))

unique(oo$h2)
unique(op$h2)

############## SCA Comparison ###########
compare_df_sca <- observed_pheno %>%
  left_join(oo, by = 'PEDIGREE_NAME') %>%
  mutate(SSGBLUP = predicted.value)
compare_df_sca <- left_join(compare_df_sca, op, by = 'PEDIGREE_NAME', suffix = c('_SSGBLUP', '_ABLUP'))

cor(compare_df_sca$TRAIT_VALUE, compare_df_sca$predicted.value_SSGBLUP, use = 'complete.obs')
cor(compare_df_sca$TRAIT_VALUE, compare_df_sca$predicted.value_ABLUP, use = 'complete.obs')
cor(compare_df_sca$predicted.value_ABLUP, compare_df_sca$predicted.value_SSGBLUP)

# reliability
mean(compare_df_sca$reliability_SSGBLUP, na.rm = T)
mean(compare_df_sca$reliability_ABLUP, na.rm = T)


############# GCA Comparison #################

# trait <- 'AVGHBRIX'
ssgblup_path <- '/mnt/Tomato/SSGBLUP'
ablup_path <- '/mnt/Tomato/ABLUP'

load(paste0(ssgblup_path, '/', trait, '_deep_fold_1', '/', 'BLUPF90_removeYearNone.RData'))

load(paste0(ablup_path, '/', trait, '_deep_fold_1', '/', 'BLUPF90_removeYearNone.RData'))

ssgblup_gca <- oo %>% 
  filter(mdl == 'P_GCA') %>% 
  select(PEDIGREE_NAME, predicted.value, reliability)

ablup_gca <- op %>% 
  filter(mdl == 'P_GCA') %>% 
  select(PEDIGREE_NAME, predicted.value, reliability)

gca_compare <- ssgblup_gca %>% 
  inner_join(ablup_gca, by = 'PEDIGREE_NAME', suffix = c('_SSGBLUP', '_ABLUP'))

cor(gca_compare$predicted.value_ABLUP, gca_compare$predicted.value_SSGBLUP)
mean(gca_compare$reliability_SSGBLUP, na.rm = T)
mean(gca_compare$reliability_ABLUP, na.rm = T)


######### Visualize GCA ###################

# SSGBLUP vs ABLUP 
parent_df <- s3read_using(read_excel, 
                          object = 'Sweetcorn/Temperate_ProcessingAndFresh/6S/Implementation/SweetCorn_ImplementationHybridList_GTStatus.xlsx',
                          bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

parent_ls <- unique(c(parent_df$P1_PEDIGREE, parent_df$P2_PEDIGREE))

inputfile_ABLUP <- '/domino/datasets/local/Testing/Check/ABLUP/'
inputfile_GBLUP <- '/domino/datasets/local/Testing/Check/Data/'

# inputfile_ABLUP <- '/domino/datasets/local/Testing/Testing/ABLUP/'
# inputfile_GBLUP <- '/domino/datasets/local/Testing/Testing/Data/'

trait <- c('SC_TF', 'EDIA', 'HSC', 'SC_HE', 'SC_SE', 'SDV', 'ESHK', 'RTLR', 'QUAL', 'SC_FG', 'S50_BE', 'LDSR_PUCCSO',
           'ELEN', 'LDSR_SETOTU','PA', 'S50', 'SC_FR', 'SC_PR', 'TILLR')

# trait <- c('SC_TF', 'EDIA', 'HSC', 'SC_HE', 'SC_SE', 'SDV', 'ESHK', 'TILLR' ,
#            'RTLR', 'QUAL', 'SC_FR','SC_FG','S50_BE', 'LDSR_SETOTU','LDSR_PUCCSO')


GCA_all <- data.frame()
for (trt in trait){
  print(trt)
  
  load(paste0(inputfile_ABLUP, trt, '_deep/BLUPF90_removeYearNone.RData'))
  load(paste0(inputfile_GBLUP, trt, '_deep/BLUPF90_removeYearNone.RData'))
  
  ABLUP_sub <- op %>% 
    filter(mdl == 'P_GCA' & PEDIGREE_NAME %in% parent_ls) %>% 
    select(PEDIGREE_NAME, predicted.value, reliability, h2)
  
  GBLUP_sub <- oo %>% 
    filter(mdl == 'P_GCA' & PEDIGREE_NAME %in% parent_ls) %>% 
    select(PEDIGREE_NAME,predicted.value, reliability, h2)
  
  GCA <- ABLUP_sub %>% 
    full_join(GBLUP_sub, by = 'PEDIGREE_NAME', suffix = c('_ABLUP', '_SSGBLUP')) %>% 
    mutate(trait = trt)
  
  if (nrow(GCA_all) == 0){
    GCA_all <- GCA
  }else{
    GCA_all <- rbind(GCA_all, GCA)
  }
  # ggplot(GCA, aes(x=predicted.value_ABLUP, y=predicted.value_SSGBLUP)) + 
  #   geom_point()
}

library(trelliscopejs)

# GCA
qplot(predicted.value_ABLUP, predicted.value_SSGBLUP, data = GCA_all) + 
  xlab('ABLUP') + 
  ylab('SSGBLUP') + 
  theme_bw() +
  facet_trelliscope(~trait, nrow = 4, ncol = 3, scales = 'free', width = 1200)


# Heritability

h2_df <- GCA_all %>% 
  select(trait, h2_ABLUP, h2_SSGBLUP) %>% 
  unique()

ggplot(h2_df, aes(x = h2_ABLUP, y = h2_SSGBLUP)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  geom_text(label = trait, hjust = 0, vjust = 0) + 
  xlab('ABLUP') + 
  ylab('SSGBLUP') + 
  theme_bw()

# Relability 

ggplot(GCA_all, aes(x = reliability_SSGBLUP)) + 
  geom_histogram() + 
  theme_bw() + 
  facet_trelliscope(~trait, nrow = 4, ncol = 3, scales = 'free', width = 1200)

# Compare reliability per trait 

reliability_df <- GCA_all %>% 
  select(PEDIGREE_NAME, reliability_ABLUP, trait) %>% 
  rename(reliability = reliability_ABLUP) %>% 
  mutate(mdl = 'ABLUP')

reliability_df_1 <- GCA_all %>% 
  select(PEDIGREE_NAME, reliability_SSGBLUP, trait) %>% 
  rename(reliability = reliability_SSGBLUP) %>% 
  mutate(mdl = 'SSGBLUP')

r2_df <- rbind(reliability_df, reliability_df_1)
r2_avg <- r2_df %>% 
  group_by(trait, mdl) %>% 
  summarize(avg_r = mean(reliability, na.rm = T))

ggplot(r2_df, aes(x = reliability, color = mdl)) + 
  geom_histogram(fill='white', position = 'dodge', bins = 20) + 
  # geom_vline(data = r2_avg, aes(xintercept = avg_r, color = mdl),
  #            linetype = 'dashed') +
  theme_bw() + 
  facet_trelliscope(~trait, nrow = 4, ncol = 3, width = 1200)



########################## Negative Reliability ####################

neg_r2 <- data.frame()

for (trt_1 in trait){
  load(paste0(inputfile_ABLUP, '/', trt_1, '_deep', '/', 'BLUPF90_removeYearNone.RData'))
  
  load(paste0(inputfile_GBLUP, '/', trt_1, '_deep', '/', 'BLUPF90_removeYearNone.RData'))
  
  ssgblup_gca <- oo %>% 
    filter(mdl == 'P_GCA') %>% 
    select(PEDIGREE_NAME, predicted.value, reliability)
  
  ablup_gca <- op %>% 
    filter(mdl == 'P_GCA') %>% 
    select(PEDIGREE_NAME, predicted.value, reliability)
  
  gca_compare <- ssgblup_gca %>% 
    inner_join(ablup_gca, by = 'PEDIGREE_NAME', suffix = c('_SSGBLUP', '_ABLUP'))
  
  # Parents that have offsprings with directly observed phenotypes 
  par_pheno <- as.character(c(Alldata$P1[Alldata$OBSRVTN_REF_CD == trt_1], Alldata$P2[Alldata$OBSRVTN_REF_CD == trt_1]))
  
  gca_compare <- gca_compare %>% 
    mutate(genotyped = if_else(PEDIGREE_NAME %in% geno$Pedigree, 1, 0)) %>% 
    mutate(directly_observed_pheno = if_else(PEDIGREE_NAME %in% par_pheno, 1, 0),
           trait = trt_1)
  
  
  temp <- gca_compare %>% 
    filter(reliability_SSGBLUP < 0)
  
  if (nrow(neg_r2) == 0){
    neg_r2 <- temp
  }else{
    neg_r2 <- rbind(neg_r2, temp)
  }
  
}

write.csv(neg_r2, file = '/mnt/neg_r2.csv', row.names = F)

