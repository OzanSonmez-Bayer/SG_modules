###################
### FALL repeated measures
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
seasons <- c('2015:08',  '2016:08', '2017:08', 
             '2018:08',  '2019:08',  '2020:08','2021:08')

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

itg_data <-dplyr::filter(itg_data, !grepl('PSL', PMC))
itg_data <-dplyr::filter(itg_data, !grepl('ASL', PMC))
itg_data <-dplyr::filter(itg_data, !grepl('BAP', PMC))

### correct field names
### filter out field names not for training set
itg_data$FIELD_NAME<-gsub("ALIR|ALL",'ALE',itg_data$FIELD_NAME)

### LEDI and LIES and other data C17 is a cucumber program but collection of data is different we have not included it in training fall
#C17<-AllData%>%filter(SET_OWNER_PROG=="CI7")
#c17_sum = C17%>%group_by(FIELD_NAME,EXPER_STAGE_REF_ID,GROWSEASON)%>%summarize(count=n())
#tst_c17<-c17_sum%>%filter(GROWSEASON%in%c("2015:08","2016:08","2017:08","2018:08","2019:08","2020:08","2021:08"))
#allC17<-unique(tst_c17$FIELD_NAME)
#C17_lies<-C17%>%filter(FIELD_NAME=="R1")

###### additional formating based on Marcel file
#itg_data$FIELD_NAME<-gsub("EHK2|K2RD|K2",'KEIJ',itg_data$FIELD_NAME)
#itg_data$FIELD_NAME<-gsub("EHA4|A4",'ALE',itg_data$FIELD_NAME)
#itg_data$FIELD_NAME<-gsub("EHH3|H3RD|H3",'HES',itg_data$FIELD_NAME)
#itg_data$FIELD_NAME<-gsub("EHR1|R1A|R1V",'LIES',itg_data$FIELD_NAME)#NO DATA
#itg_data$FIELD_NAME<-gsub("EHB3",'BAK',itg_data$FIELD_NAME)
#itg_data$FIELD_NAME<-gsub("EHB4",'BOX',itg_data$FIELD_NAME)
#itg_data$FIELD_NAME<-gsub("EHH1|H1|H1RD",'HAK',itg_data$FIELD_NAME)
#itg_data$FIELD_NAME<-gsub("EHH2|H2",'HELM',itg_data$FIELD_NAME)
#itg_data$FIELD_NAME<-gsub("EHL1|L1",'LEDI',itg_data$FIELD_NAME)

field_names<-c("ALE","BAK","BOX","HAK","HELM","HES","KEIJ","LEDI","LIES","V29","V30","V31","V41")
itg_data = itg_data %>% filter(FIELD_NAME%in%field_names)
all_fall<-itg_data

### unit conversion   
all_fall <- conv_uom(all_fall) 
all_fall<-all_fall%>%filter(!VALUE_TYPE == "Text")
all_fall<-all_fall%>%filter(!VALUE_TYPE == "Character list")

### subset by traits
pheno=all_fall
TRAITS<-c("EXTCO","FRLGT","AFW_C","FRNMK","SHAPE","NETWT","FRNMK","NHPLT","NPA" )
pheno<-pheno%>%filter(OBSRVTN_REF_CD%in%TRAITS)

### Outliers after visual plots evaluation *** this should be changed to use distribution not per test set but per raw dataset
dat=pheno %>%
  dplyr::group_by(OBSRVTN_REF_CD)%>%
  dplyr::summarise(min=min(na.omit(TRAIT_VALUE)),max=max(na.omit(TRAIT_VALUE)),median=median(na.omit(TRAIT_VALUE)))%>%
  data.frame()


##outlier FRLGT 
frl=pheno%>%filter(OBSRVTN_REF_CD=="FRLGT" & TRAIT_VALUE >60)
pheno[rownames(pheno) %in% rownames(frl),]$TRAIT_VALUE=dat1[which(dat$OBSRVTN_REF_CD=="FRLGT"),"median"]

##outlier FRNMK 
#frn=pheno%>%filter(OBSRVTN_REF_CD=="FRNMK" & TRAIT_VALUE >40)
#pheno[rownames(pheno) %in% rownames(frl),]$TRAIT_VALUE=dat[which(dat$OBSRVTN_REF_CD=="FRLGT"),"median"]

##outlier NETWT
#net=pheno%>%filter(OBSRVTN_REF_CD=="NETWT" & TRAIT_VALUE >15)
#pheno[rownames(pheno) %in% rownames(frl),]$TRAIT_VALUE=dat[which(dat$OBSRVTN_REF_CD=="FRLGT"),"median"]


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
save(pheno,file="Cln_cucumber_fall12102021_new.Rdata")

