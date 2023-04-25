###################
### Summer dataset
###################
source("/repos/calc_trt_git.R")# change this name to your repo
source("/repos/sec_vault_git.R")# change this name to your repo
load_env()
#### parameters ####
#@Crop can be modified with any name
#@seasons 05 summer, 08 fall,01 spring
bucket <- 'veg-apd-sdi-predictiveanalytcs-prod-pheno-data'
path <- '/H2H_Data/'
Crop <- 'Cucumber'
prog <- 'Y3'
seasons <- c('2015:05',  '2016:05', '2017:05', 
             '2018:05',  '2019:05',  '2020:05','2021:05')
years <- c('2014','2015','2016','2017','2018','2019','2020','2021')

## Current training set modify accordingly as coding of lines can change pedigree names and us updated weekly
## in the system also used in ODB
#### Pull Phenotypic data ####
AllData <- data.frame()
for (selectedYear in years){
  infile <- paste0(path,Crop,"/",selectedYear,"/AllYearData.RData")
  s3load(object = infile, bucket = "veg-apd-sdi-predictiveanalytcs-prod-pheno-data")
  AllData <- rbind(AllData,YearPhenos)
}

#### FILTERING ####
# this is done in consultation with the breeder and assistant (Carin Veen)
AllData<-  AllData%>%filter(SET_OWNER_PROG %in% c("Y3","NR5") & GROWSEASON%in%seasons)
itg_data <- AllData %>% 
  filter(GROWSEASON %in% seasons & EXPER_STAGE_REF_ID%in%c("P1","P2","P3","P4","S1"))
itg_data <-dplyr::filter(itg_data, !grepl('BAP', PMC))
itg_data <-dplyr::filter(itg_data, !grepl('HIGH WIRE', PMC))
itg_data <-dplyr::filter(itg_data, !grepl('HW', TEST_SET_NAME))
itg_data = itg_data%>%filter(!TEST_SET_NAME=="21BAPP1A51")

### correct field names
### filter out field names not for training set
itg_data$FIELD_NAME<-gsub("ALIR|ALL",'ALE',itg_data$FIELD_NAME)
itg_data = itg_data%>%filter(!FIELD_NAME%in%c("CGM","VGE","BAP", "V30","CLSC","V42"))

### check on repetition by test set and verify is within expection per season
dat_rep = itg_data%>%group_by(OBSRVTN_REF_CD,GROWSEASON,TEST_SET_NAME)%>%dplyr::summarise(count=max(REPETITION),val=unique(PMC))
all_summer<-itg_data

### unit conversion      
all_summer <- conv_uom(all_summer) 
all_summer<-all_summer%>%filter(!VALUE_TYPE == "Text")
all_summer<-all_summer%>%filter(!VALUE_TYPE == "Character list")
all_summer$TRAIT_VALUE = as.numeric(all_summer$TRAIT_VALUE)

### remove P3 distribution of traits is different
pheno=all_summer%>% filter(EXPER_STAGE_REF_ID%in%c("P1","P2","S1"))
 
### Outliers after visual plots evaluation *** this should be changed to use distribution not per test set but per raw dataset
##outlier NETWT
dat=pheno %>%dplyr::group_by(OBSRVTN_REF_CD)%>%dplyr::summarise(min=min(na.omit(TRAIT_VALUE)),max=max(na.omit(TRAIT_VALUE)),median=median(na.omit(TRAIT_VALUE)))%>%data.frame()
net=pheno%>%filter(OBSRVTN_REF_CD=="NETWT" & TRAIT_VALUE >15)
pheno[rownames(pheno) %in% rownames(net),]$TRAIT_VALUE=dat[which(dat$OBSRVTN_REF_CD=="NETWT"),"median"]
 
### subset by traits

TRAITS<-c("EXTCO","FRLGT","AFW_C","FRNMK","SHAPE","NETWT","FRNMK","NHPLT","NPA" )
pheno<-pheno%>%filter(OBSRVTN_REF_CD%in%TRAITS)

#### CALCULATED TRAITS ####
set<-unique(pheno$TEST_SET_NAME)
pheno_all=data.frame()
for(i in 1:length(set)){
  pheno1 = pheno%>%filter(TEST_SET_NAME==set[i])
  tst=calc_trt(pheno=pheno1,Crop,prog)
  pheno_all=rbind(pheno_all,tst)
}


#### REPEATED MEASURES ####
vals=repps(pheno_all)
avg_all = Avg_row_head(vals)
sum_all = Sum_rep(vals)
pheno1<-avg_all[["Pheno"]]
pheno2<-sum_all[["Pheno"]]
pheno3 = pheno2%>%filter(!OBSRVTN_REF_CD%in%pheno1$OBSRVTN_REF_CD)
pheno=rbind(pheno3,pheno1)

#### save results I include date as data in the system updates frequently
save(pheno,file="Cln_cucumber_summer18022022_new.Rdata")

