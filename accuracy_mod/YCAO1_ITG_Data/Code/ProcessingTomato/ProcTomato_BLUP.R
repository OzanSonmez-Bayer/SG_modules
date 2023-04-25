setwd('C:/Users/ycao1/OneDrive - Monsanto/Migrated from My PC/Desktop/VEG_ITG')

library(asreml)
library(tidyverse)

dat <- read_csv('PROCESSING_TOMATO_9Z_PLOT_OBS.csv')

# Clean data
# Use AVGHBRIX as an example 
dat_clean <- dat %>% 
  select('Cross Name', 'Field', 'Origin', 'Origin (Female Parent)', 'Origin (Male Parent)', 'Pedigree', 
         'Pedigree (Female Parent)', 'Pedigree (Male Parent)', 'Rep#', 'Season','Set', 
         'PO.AVGHBRIX.NOT APPLICABLE.NOT APPLICABLE.01')
names(dat_clean) <- c('CrossName', 'Field', 'Origin', 'Origin_female', 'Origin_male', 'Pedigree', 'Pedigree_female', 
                      'Pedigree_male', 'Rep', 'Season', 'Set', 'AVGHBRIX')

# Season: 2018:04, 2018:01, 2018
# NA in Rep

# User season 2018:04 and 2018:01
hbrix<- dat_clean %>% 
  filter(Season %in%  c( "'2018:04", "'2018:01")) %>% 
  filter(Field != 'ALEX') %>% 
  filter(Set != 'ALEX')

# Pedigree based model
# NAs in Pedigree_female and Pedigree_male
pd <- hbrix %>% 
  select(Pedigree, Pedigree_female, Pedigree_male) %>% 
  unique() %>% 
  mutate(Pedigree_female = if_else(is.na(Pedigree_female), Pedigree, Pedigree_female)) %>% 
  mutate(Pedigree_male = if_else(is.na(Pedigree_male), Pedigree, Pedigree_male))

# Pedigree Relationship Matrix
names(pd) <- c('Pedigree', 'Female', 'Male')
pd.ainv <- asreml.Ainverse(pd) # Erros, couldn't get A matrix. Use GCA model

# GCA model
fm_hbrix <- asreml(AVGHBRIX~Season + Field + Rep + Set, random = ~Pedigree_female + and(Pedigree_male),
                   data = hbrix, equate.levels = c('Pedigree_female', 'Pedigree_male'), na.method.X = 'include')

# BLUP
blup_hbrix <- data.frame(round(coef(fm_hbrix)$random,4), round(fm_hbrix$vcoeff$random,4))
names(blup_hbrix) <- c('BLUP', 'BLUP.se')

# Results don't look correct 

# AVJB
dat_avjb <- dat %>% 
  select('Cross Name', 'Field', 'Origin', 'Origin (Female Parent)', 'Origin (Male Parent)', 'Pedigree', 
         'Pedigree (Female Parent)', 'Pedigree (Male Parent)', 'Rep#', 'Season','Set', 
         'PO.AVJB.NOT APPLICABLE.NOT APPLICABLE.01')
names(dat_avjb) <- c('CrossName', 'Field', 'Origin', 'Origin_female', 'Origin_male', 'Pedigree', 'Pedigree_female', 
                      'Pedigree_male', 'Rep', 'Season', 'Set', 'AVJB')

avjb<- dat_avjb %>% 
  filter(Season %in%  c( "'2018:04", "'2018:01")) %>% 
  filter(Field != 'ALEX') %>% 
  filter(Set %in% c('2NDEARLYI3', 'AGSD_EARLY','AGSD_F3', 'AGSD_MAIN', 'AO_EXTRA', 'AO_THICK17', '17SE_I3', '17SE_I3PEN',
                    '18SE1', '18SE2', '18SE3'))

fm_avjb <- asreml(AVJB~Season + Field + Rep + Set, random = ~Pedigree_female + and(Pedigree_male),
                   data = avjb, equate.levels = c('Pedigree_female', 'Pedigree_male'), na.method.X = 'include')

# BLUP
blup_avjb <- data.frame(round(coef(fm_avjb)$random,4), round(fm_avjb$vcoeff$random,4))
names(blup_avjb) <- c('BLUP', 'BLUP.se')
blup_avjb$Pedigree <- gsub('Pedigree_female_', '', rownames(blup_avjb))


# Function to calculate GCA
# Exclue season 2018 and field ALEX, Set ALEX

blupFun <- function(data){
  bdat <- data %>%
    filter(Season %in%  c( "'2018:04", "'2018:01")) %>% 
    filter(Field != 'ALEX') %>% 
    filter(Set != 'ALEX')
  
  bdat <- bdat %>%
    select(Pedigree, Pedigree_female, Pedigree_male) %>% 
    unique() %>% 
    mutate(Pedigree_female = if_else(is.na(Pedigree_female), Pedigree, Pedigree_female)) %>% 
    mutate(Pedigree_male = if_else(is.na(Pedigree_male), Pedigree, Pedigree_male))
  
  fm_mdl <- asreml(Trait_value~Season + Field + Rep + Set, random = ~Pedigree_female + and(Pedigree_male),
                   data = bdat, equate.levels = c('Pedigree_female', 'Pedigree_male'), na.method.X = 'include')
  blup_df <- data.frame(round(coef(fm_mdl)$random,4), round(fm_mdl$vcoeff$random,4))
  names(blup_df) <- c('BLUP', 'BLUP.se')
  
  blup_df$Pedigree <- gsub('Pedigree_female_', '', rownames(blup_df), fixed = T)
  
  # pedigree df
  p_df <- table(c(Pedigree_female, Pedigree_male))
  
  blup_df$Num_Hybrid <- unname(p_df[blup_df$Pedigree])
  return(blup_df)
}

# AVGHBRIX
dat_avghbrix <- dat %>% 
  select('Cross Name', 'Field', 'Origin', 'Origin (Female Parent)', 'Origin (Male Parent)', 'Pedigree', 
         'Pedigree (Female Parent)', 'Pedigree (Male Parent)', 'Rep#', 'Season','Set', 
         'PO.AVGHBRIX.NOT APPLICABLE.NOT APPLICABLE.01')
names(dat_avghbrix) <- c('CrossName', 'Field', 'Origin', 'Origin_female', 'Origin_male', 'Pedigree', 'Pedigree_female', 
                      'Pedigree_male', 'Rep', 'Season', 'Set', 'Trait_value')
blup_hbrix1 <- blupFun(dat_avghbrix)
