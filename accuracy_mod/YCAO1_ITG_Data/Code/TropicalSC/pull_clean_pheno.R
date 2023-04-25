library(aws.s3)
library(httr)
library(tidyverse)
install.packages("azurequest", repos=c('https://cran.science-at-scale.io', 'https://cran.rstudio.com'))
library(azurequest)
library(asreml)
library(asremlPlus)
install.packages('packages/doBy_4.6-2.tar.gz', type = 'source', repos = NULL)
library(doBy)
library(measurements)
library(outliers)



# source('VaultCreds.R')
# source('ValutCreds_ancestry.R')

source('VaultCredentials.R')
VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")



seasons <- c('2018:07', '2018:09', '2019:03', '2019:07', '2019:09','2020:03', '2020:07', '2020:09', '2021:03')

Crop <- 'Corn'

bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
# years <- c('2018', '2019', '2020')

years <- c('2018', '2019', '2020', '2021')
AllData <- data.frame()
for (selectedYear in years){
  infile <- paste0(path,Crop,"/",selectedYear,"/AllYearData.RData")
  s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")
  if (nrow(AllData) !=0){
    YearPhenos <- YearPhenos %>% 
      select(colnames(AllData))
    AllData <- rbind(AllData,YearPhenos) 
  }else{
    AllData <- YearPhenos
  }
}

AllData_sub <- AllData %>% 
  filter(GROWSEASON %in% seasons) 


# %>% 
#   filter(SET_OWNER_PROG == 'AV') 


# %>% 
#   mutate(OBSRVTN_REF_CD = paste0(OBSRVTN_REF_CD,'_', DESCRIPTOR))

AllData_sub <- AllData_sub[grep('TSY-M-PROCESSING & FRESH', as.character(AllData_sub$PMC)), ]

AllData_sub_1 <- AllData_sub %>% 
  filter(DESCRIPTOR %in% c('SPIRKU', 'PHSPMA', 'SETOTU')) %>% 
  mutate(OBSRVTN_REF_CD = paste0(OBSRVTN_REF_CD, '_', DESCRIPTOR))


AllData_sub_2 <- AllData_sub %>% 
  filter(!DESCRIPTOR %in% c('SPIRKU', 'PHSPMA', 'SETOTU'))


AllData_sc <- rbind(AllData_sub_1, AllData_sub_2)

unique(AllData_sc$OBSRVTN_REF_CD)
unique(AllData_sc$GROWSEASON)
unique(AllData_sc$PMC)

AllData_sc <- AllData_sc %>% 
  filter(OBSRVTN_REF_CD %in% c('EDIA', 'EHT', 'ELEN', 'ERSC', 'EVG', 'HSC', 'KRN', 
                               'LDSR_SPIRKU', 'LDSR_PHSPMA', 'LDSR_SETOTU', 'P50D', 'PCTRC', 
                               'PHT', 'QUAL', 'RTLP', 'S50D', 'SC_FR', 'SC_HE', 'SC_TF', 'SCCSA', 'SCTPA'))


##################### Outlier Analysis ########################

outL = function(X){
  if (length(X)==0){out=NULL}else{
    pval = grubbs.test(X)$p.value
    if (pval>0.01 | is.na(pval)){ out=X }else{
      out = rm.outlier(X, fill=TRUE,median=TRUE)
    }
  }
  list(out = out, pval = pval)
}

# for a specific trait and set combo run the outlier algorithm and replace it
# by the median and spit out the new data set
rm.out = function(pheno,PRG,tst='Grub'){
  #Data = subset(pheno, SET_OWNER_PROG==PRG)
  trait = unique(pheno$OBSRVTN_REF_CD)
  test = unique(pheno$TEST_SET_NAME)
  pheno$TRAIT_VALUE=as.numeric(as.character(pheno$TRAIT_VALUE))
  #data.new = pheno
  Tab = data.frame(Row=c(),Cross = c(), SetName = c(),PlotBid=c(),PlotNumber = c(), Trait =  c(), Season =  c(),
                   Pvalue =  c(), TraitVal = c(), mean =  c(), 
                   median =  c(), SD =  c(),Samples=c(), Series=c())
  
  for (i in 1:length(trait)){
    for (j in 1: length(test)){
      dat = subset(pheno, TEST_SET_NAME==test[j] & OBSRVTN_REF_CD==trait[i])
      if (dim(dat)[1]!=0){
        ###Grubs test
        if (tst=='Grub'){
          tesT = outL(as.numeric(dat$TRAIT_VALUE))
          c = tesT$out
          p = tesT$pval
        }
        ###MAD
        if (tst=='MAD'){
          tesT = mad_based_outlier(dat,'TRAIT_VALUE',thresh=12)
          c=tesT[[2]]
          p = tesT[[1]]
        }
        if (sum(is.na(c))==length(c)) c=as.numeric(dat$TRAIT_VALUE)
        out.w = which(c!=dat$TRAIT_VALUE)
        m = length(out.w)
        data.out = dat[out.w, ]
        Mean = round(mean(dat$TRAIT_VALUE),2)
        Median = round(median(dat$TRAIT_VALUE),2)
        SD = round(sd(dat$TRAIT_VALUE),2)
        Sample = paste(dat$TRAIT_VALUE, collapse = ",")
        Series=paste(hist(dat$TRAIT_VALUE,plot=F)$counts,collapse=",")
        # PlotNumber <- paste(unique(data.out$PLOT_NUMBER), collapse = ',')
        V = data.frame(Row=rownames(data.out),Cross=data.out$CROSS_NAME, 
                       SetName = rep(test[j], m), PlotBid=data.out$PLOT_BID,
                       PlotNumber =  rep(paste(unique(data.out$PLOT_NUMBER), collapse = ','),m),
                       Trait = rep(trait[i], m), 
                       Season = data.out$GROWSEASON,Pvalue = rep(round(p,5), m), 
                       TraitVal = round(data.out$TRAIT_VALUE,2), mean = rep(Mean, m), 
                       median = rep(Median, m), SD = rep(SD,m),Sample = rep(Sample,m),Series=rep(Series,m))
        Tab = rbind(Tab, V)
        #dat$TRAIT_VALUE = c
        #data.new[data.new$TEST_SET_NAME==test[j] & data.new$OBSRVTN_REF_CD==trait[i], ] = dat
      }
    }
  }
  Tab
}



outliers <- rm.out(AllData_sc, PRG = 'AV')
outliers <- outliers %>% 
  mutate(ID = paste0(Row, '_', Cross, '_', SetName, '_', PlotBid, '_', PlotNumber, '_', Trait, '_', Season))

AllData_sc_cln <- AllData_sc %>% 
  mutate(ID = paste0(V1,'_', CROSS_NAME, '_', TEST_SET_NAME,  '_', PLOT_BID, '_', PLOT_NUMBER, '_', OBSRVTN_REF_CD, '_', GROWSEASON)) %>% 
  filter(!ID %in% outliers$ID)
AllData_sc_cln <- AllData_sc_cln%>% 
  select(-ID)


s3saveRDS(AllData_sc_cln, object = 'ycao1/AoB/SweetCorn/Tropical/rawdata/pheno_2018_2021_clean_08112021.rds',
          bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
