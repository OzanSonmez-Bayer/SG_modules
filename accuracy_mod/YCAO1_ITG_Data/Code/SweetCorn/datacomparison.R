# Compare data pulled based on PMC and Sets

install.packages('aws.s3')
library(aws.s3)
library(dplyr)
library(jsonlite)
library(plyr)
library(readxl)
library(asreml)
library(asremlPlus)
library(doBy)
library(reshape)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggrepel)

source('/repos/ITG_Data/Code/Shared/misFunctions.R')

AWS.Authorization('ycao1')

path <- 'ycao1/ITG_data/SweetCorn/Data/'
bucket <- 'genome-analytics-perm-space'

fresh_data <- data.frame()
proc_data <- data.frame()

file_list <- get_bucket(bucket = bucket, prefix = paste(path, 'Raw/', sep = ""))
file_list <- file_list[grepl('_Raw', file_list)]

for (i in 1:length(file_list)){
  print(file_list[[i]]$Key)
  infile <- paste0('s3://','genome-analytics-perm-space/', file_list[[i]]$Key)
  csvfile <- get_object(infile)
  csvfile_1 <- rawToChar(csvfile)
  con <- textConnection(csvfile_1)
  temp_dat <- read.csv(con)
  close(con)
  if (grepl('_Fresh_', file_list[[i]]$Key)){
    fresh_data<- plyr::rbind.fill(fresh_data, temp_dat)
  }else{
    proc_data <- plyr::rbind.fill(proc_data, temp_dat)
  }
}


fresh_pmc <- s3read_using(read.csv, 
                          object = "s3://genome-analytics-perm-space/ycao1/Horizon/Corn/Corn_FRESH.csv")
proc_pmc <- s3read_using(read.csv, 
                         object = "s3://genome-analytics-perm-space/ycao1/Horizon/Corn/Corn_PROCESSING_SEY-M-PROCESSING.csv")

# Subset fresh_pmc with for year 2015:05, 2016:04, 2017:04, and 2018:04

fresh_pmc_subset <- fresh_pmc %>% 
  filter(bmrd.grow_year %in% c(2015, 2016, 2017, 2018) & bmrd.grow_month == 4 & bmrd.grower == '6S')

proc_pmc_subset <- proc_pmc %>% 
  filter(bmrd.grow_year %in% c(2015, 2016, 2017, 2018) & bmrd.grow_month == 4 & bmrd.grower == '6S')


fresh_set <- fresh_data %>% 
  filter(GROWSEASON %in% c('2015:04', '2016:04', '2017:04', '2018:04'))


proc_set <- proc_data %>% 
  filter(GROWSEASON %in% c('2015:04', '2016:04', '2017:04', '2018:04'))


# Compare BLUPs (PBLUP)
# Processing corn as an example

proc_pmc_subset <- proc_pmc %>% 
  filter(bmrd.grow_year %in% c(2017, 2018, 2019) & bmrd.grow_month == 4)

proc_pmc_sc123 <- proc_pmc_subset %>% 
  filter(bmrd.stage %in% c('S1', 'S2', 'S3'))

length(unique(proc_pmc_sc123$bmrd.set_name))

# BLUPs with mixed sets

source('/repos/ITG_Data/Code/Shared/BLUPFunctions.R')

sets <- read_excel('/repos/ITG_Data/Code/SweetCorn/dataPUll_PMC.xlsx', sheet = 'SC123_PROC')
colnames(sets) <- c('TEST_SET_NAME', 'Decision')

blup_w_mixed_sets <- pblup_fun(proc_pmc_sc123, crop = 'Corn', pmc = 'processing', year = NA)
pa_w_mixed_sets <- pblup_fun(proc_pmc_sc123[proc_pmc_sc123$bord.obsrvtn_ref_cd == 'PA',], 
                             crop = 'Corn', pmc = 'processing', year = NA)

proc_pmc_sc123_1 <- proc_pmc_sc123 %>% 
  filter(.$bmrd.set_name %in% sets$TEST_SET_NAME[is.na(sets$Decision)])


blup_wo_mixed_sets <- pblup_fun(proc_pmc_sc123_1, crop = 'Corn', pmc = 'processing', year = NA)

pa_no_mixed_sets <- pblup_fun(proc_pmc_sc123_1[proc_pmc_sc123_1$bord.obsrvtn_ref_cd == 'PA',], 
                              crop = 'Corn', pmc = 'processing', year = NA)

write.csv(blup_w_mixed_sets, file = '/mnt/BLUP_PROC_MixedSets.csv')
write.csv(blup_wo_mixed_sets, file = '/mnt/BLUP_PROC_no_MixedSets.csv')


BLUP_w_mixed_sets <- read.csv('/repos/ITG_Data/Code/SweetCorn/BLUP_PROC_MixedSets.csv')
blup_wo_mixed_sets <- read.csv('/repos/ITG_Data/Code/SweetCorn/BLUP_PROC_no_MixedSets.csv')

BLUP_w_mixed_sets_1 <- BLUP_w_mixed_sets[-grep('^PA', BLUP_w_mixed_sets$traitRefId),-1]
BLUP_w_mixed_sets <- rbind(BLUP_w_mixed_sets_1, pa_w_mixed_sets)

blup_wo_mixed_sets_1 <- blup_wo_mixed_sets[-grep('^PA', blup_wo_mixed_sets$traitRefId),-1]
blup_wo_mixed_sets <- rbind(blup_wo_mixed_sets_1, pa_no_mixed_sets)

# Cast
blup_w_sca <- BLUP_w_mixed_sets[grep('_PSCA_', BLUP_w_mixed_sets$traitRefId),]

blup_w_sca <- blup_w_sca %>% 
  dplyr::select(lineName, traitRefId, value) %>% 
  spread(traitRefId, value)

blup_w_gca <- BLUP_w_mixed_sets[grep('_PGCA_', BLUP_w_mixed_sets$traitRefId),]

blup_w_gca <- blup_w_gca %>% 
  dplyr::select(lineName, traitRefId, value) %>% 
  spread(traitRefId, value)

blup_wo_sca <- blup_wo_mixed_sets[grep('_PSCA_', blup_wo_mixed_sets$traitRefId),] 
blup_wo_sca <- blup_wo_sca %>% 
  dplyr::select(lineName, traitRefId, value) %>% 
  spread(traitRefId, value)

blup_wo_gca <- blup_wo_mixed_sets[grep('_PGCA_', blup_wo_mixed_sets$traitRefId),] 
blup_wo_gca <- blup_wo_gca %>% 
  dplyr::select(lineName, traitRefId, value) %>% 
  spread(traitRefId, value)



# blup_w_sca <- blup_w[, grep('_PSCA_', colnames(blup_w))]
# blup_w_sca$lineName <- blup_w$lineName
# 
# miss_num_w_sca <- apply(blup_w_gca, MARGIN = 1, FUN = function(x) {sum(is.na(x))})
# 
# blup_w_gca <- blup_w[, grep('_PGCA_', colnames(blup_w))]
# blup_w_gca$lineName <- blup_w$lineName
# 
# blup_wo_sca <- blup_wo[, grep('_PSCA_', colnames(blup_wo))]
# blup_wo_sca$lineName <- blup_wo$lineName
# 
# blup_wo_gca <- blup_wo[, grep('_PGCA_', colnames(blup_wo))]
# blup_wo_gca$lineName <- blup_wo$lineName

# Compare blups

sca_compare <- blup_wo_sca %>% 
  dplyr::inner_join(blup_w_sca, by = 'lineName', suffixes = c("_wo","_w"))


gca_compare <- blup_wo_gca %>% 
  dplyr::inner_join(blup_w_gca, by = 'lineName', suffixes = c("_wo","_w"))

# Reliability
trtlist <- c("EDIA" , "ELEN", "HSC","PA", "PCTRC","RTLR" , "SC_HE" ,"SC_TF" ,"SCCSA", "SCTPA", "SCUNH")

dat_comp <- sca_compare
# Average reliability

dat_comp[,-1] <- apply(dat_comp[,-1], 2, as.numeric )
reliability <- colMeans(dat_comp[, grep('_RELIABILITY', colnames(dat_comp))], na.rm = T)

reliability_df <- data.frame('Trait' = names(reliability)[1:11], 
                              'With_Mixed_Sets' = unname(reliability)[1:11],
                              'Without_Mixed_Sets' = unname(reliability)[12:22])
reliability_df$Trait <- unlist(lapply(as.character(reliability_df$Trait),
                                       FUN = function(x) {unlist(strsplit(x, split = '_'))[1]}))
reliability_df$Trait[c(7,8)] <- c('SC_HE', 'SC_TF')

reliability_df

ggplot(reliability_df, aes(x = With_Mixed_Sets, y = Without_Mixed_Sets, label = Trait)) + 
  geom_point(color = 'blue') + 
  geom_text_repel(aes(label = Trait), size = 3.5) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_minimal() + 
  labs(title = 'Reliability of GCA Comparison: ')


for(trt in trtlist){
  temp_x <- as.numeric(dat_comp[, paste0(trt, '_PGCA_RELIABILITY.x')])
  temp_y <- as.numeric(dat_comp[, paste0(trt, '_PGCA_RELIABILITY.y')])
  
  temp <- data.frame('Reliability' = c(temp_x, temp_y),
                     'Data' =  rep(c('With Mixed Sets', 'Without Mixed Sets'), c(length(temp_x), length(temp_y))))
  print(trt)
  print(cor(temp_x, temp_y, method = 'spearman', use = 'complete.obs'))
                                   
  p <- ggplot(temp, aes(x = Data, y = Reliability, fill = Data)) + 
    geom_violin() + 
    stat_summary(fun.y = mean, geom = 'point', color = "blue", size = 5) +
    labs(title = paste("Reliability: ", trt)) + 
    xlab("Data")  
  p
}

# cor(dat_comp$EDIA_PSCA_RELIABILITY.x, dat_comp$EDIA_PSCA_RELIABILITY.y, use = "complete.obs")


# Heritability
heritability <- colMeans(dat_comp[, grep('_HERITABILITY', colnames(dat_comp))], na.rm = T)
  
heritability_df <- data.frame('Trait' = names(heritability)[1:11], 
                              'With_Mixed_Sets' = unname(heritability)[1:11],
                              'Without_Mixed_Sets' = unname(heritability)[12:22])
heritability_df$Trait <- unlist(lapply(as.character(heritability_df$Trait),
                                       FUN = function(x) {unlist(strsplit(x, split = '_'))[1]}))
heritability_df$Trait[c(7,8)] <- c('SC_HE', 'SC_TF')
heritability_df$With_Mixed_Sets[heritability_df$Trait == 'PA'] <- 0.0798
heritability_df$Without_Mixed_Sets[heritability_df$Trait == 'PA'] <- 0.0668

ggplot(heritability_df, aes(x = With_Mixed_Sets, y = Without_Mixed_Sets, label = Trait)) + 
  geom_point(color = 'blue') + 
  geom_text_repel(aes(label = Trait), size = 3.5) + 
  geom_abline(intercept = 0, slope = 1) + 
  theme_minimal() + 
  labs(title = 'Heritability Comparison')

# Compare BLUPs

for(trt in trtlist){
  temp_x <-dat_comp[, paste0(trt, '_PSCA_PRED.x')]
  temp_y <-dat_comp[, paste0(trt, '_PSCA_PRED.y')]
  
  temp <- data.frame('PRED_With_Mixed_Sets' = temp_x,
                     'PRED_without_Mixed_Sets' = temp_y)
  rank_cor <- cor(temp_x, temp_y, use = 'complete.obs', method = 'spearman')
  print(rank_cor)
  
  p <- ggplot(temp, aes(x = PRED_With_Mixed_Sets, y = PRED_without_Mixed_Sets)) + 
    geom_point(color = 'blue') + 
    geom_abline(intercept = 0, slope = 1) + 
    theme_minimal() + 
    labs(title = paste('Predicted Value Comparison: ', trt)) +
    annotate("text", x = min(temp_x, na.rm = T), y = max(temp_y, na.rm = T),
             label = "paste(italic(R) ^ 2, \" = 1\")", parse = TRUE)
  p
}

# Rank correlation for top 10% BLUP
dat_store <- data.frame()
miss_ind <- list()
for(trt in trtlist){
  temp_x <- as.numeric(dat_comp[, paste0(trt, '_PGCA_PRED.x')])
  temp_y <- as.numeric(dat_comp[, paste0(trt, '_PGCA_PRED.y')])
  temp_x_r <- as.numeric(dat_comp[, paste0(trt, '_PGCA_RELIABILITY.x')])
  temp_y_r <- as.numeric(dat_comp[, paste0(trt, '_PGCA_RELIABILITY.y')])
  temp_x_n <- as.numeric(dat_comp[, paste0(trt, '_PGCA_X.x')])
  temp_y_n <- as.numeric(dat_comp[, paste0(trt, '_PGCA_X.y')])

  temp <- data.frame('lineName' = as.character(dat_comp$lineName),
                     'PRED_With_Mixed_Sets' = temp_x,
                     'PRED_without_Mixed_Sets' = temp_y,
                     'RELIABILITY_With_Mixed_Sets' = temp_x_r,
                     'RELIABILITY_without_Mixed_Sets' = temp_y_r,
                     'N_With_Mixed_Sets' = temp_x_n,
                     'N_without_Mixed_Sets' = temp_y_n,
                     'Trait' = trt)

  # temp <- temp[order(temp$PRED_without_Mixed_Sets, decreasing = T),]
  temp <- temp %>%
    filter(!is.na(PRED_With_Mixed_Sets)) %>%
    filter(!is.na(PRED_without_Mixed_Sets))
  
  temp <- temp %>% 
    arrange(desc(PRED_With_Mixed_Sets)) %>% 
    mutate(rank_with_mixed_sets = 1:nrow(.)) 
  
  temp_sub_1 <- temp[1:ceiling(nrow(temp)*0.15),]
  
  temp_1 <- temp%>% 
    arrange(desc(PRED_without_Mixed_Sets)) %>% 
    mutate(rank_without_mixed_sets = 1:nrow(.)) 
  
  temp_sub <- temp_1[1:ceiling(nrow(temp)*0.15),]
  
  print(temp_sub_1$lineName[which(!temp_sub_1$lineName %in% intersect(temp_sub_1$lineName, temp_sub$lineName))])
  
  miss_temp <- as.character(temp_sub_1$lineName[which(!temp_sub_1$lineName %in% intersect(temp_sub_1$lineName, temp_sub$lineName))])
  miss_ind <- c(miss_ind, list(miss_temp))
  
  dat_store <- rbind(dat_store, temp_sub)
  rank_cor <- cor(temp_sub$PRED_With_Mixed_Sets, temp_sub$PRED_without_Mixed_Sets, 
                  use = 'complete.obs', method = 'spearman')
  print(trt)
  print(rank_cor)
  
  # p <- ggplot(temp_sub, aes(x = PRED_With_Mixed_Sets, y = PRED_without_Mixed_Sets)) + 
  #   geom_point(color = 'blue') + 
  #   geom_abline(intercept = 0, slope = 1) + 
  #   theme_minimal() + 
  #   labs(title = paste('Top 10% of Predicted Value Comparison: ', trt)) +
  #   annotate("text", x = min(temp_sub$PRED_With_Mixed_Sets, na.rm = T), 
  #            y = max(temp_sub$PRED_without_Mixed_Sets, na.rm = T),
  #            label = "paste(italic(R) ^ 2, \" = 0.82\")", parse = TRUE)
  # p
}
# names(miss_ind) <- trtlist

df3 <- purrr::map(miss_ind, data.frame) %>%
  set_names(trtlist) %>% 
  bind_rows(.id = "trt")

tst_dat <- blup_w_gca
tst <- tst_dat[, grep('_PGCA_PRED$', colnames(tst_dat))]
tst$lineName <- tst_dat$lineName
apply(tst, MARGIN = 2, function(x){length(unique(tst$lineName[!is.na(x)]))})  
