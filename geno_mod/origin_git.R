##################
# generalized code
##################
source("/repos/ELBFA_GS_modules/geno_mod/func_sql.R")
source("/repos/ELBFA_GS_modules/geno_mod/sec_vault_git.R")
source("/repos/ELBFA_GS_modules/geno_mod/util_funcs_geno.R")

load_env()
season="summer"
year="2020"

#This code takes as input
#genotype inbred files in HD and processed 0,1,2 format, at the moment I input every inbred file available 
#a reference PCM1 file which should contain the hybrids to predict, usually current season DH in combination with testers
#the reference PCM1 file can contain more than one year and/or season that is why there is a filtering step, otherwise use just the current season PCM1 file
#


#### read file ####
#@ provided by breeder with hybrid pedigree, female invbid,male tester invbid,female and male pedigree
#dhSP<-read_excel("/repos/ELBFA_GS_modules/geno_mod/Files_PCM/List of testers of PCM1 by season by year_InBID.xlsx",sheet = 2)
dhSP<- s3read_using(FUN = read_excel, object = 's3://veg-apd-sdi-predictiveanalytcs-prod-workspace/ycao1/ITG/Cucumber/Files_PCM/Files_PCM/List of testers of PCM1 by season by year_InBID.xlsx', sheet = 2)

dhSP<-dhSP%>%dplyr::filter(Season==season & Year==year)

#### Parents pedigrees
dhSP1=data.frame(Inv_BID=c(dhSP$Female_InvBID,dhSP$`Male tester_InvBID`),Pedigree=c(dhSP$Female_pedigree,dhSP$`Male tester_pedigree`))

#### Hybrids Parents Pairs
dhSP2<-dhSP%>%select(pedigree,Female_InvBID,`Male tester_InvBID`)

#### Lists for sql
# search by invBID
invbid=data.frame(id=unique(c(dhSP$Female_InvBID,dhSP$`Male tester_InvBID`)))# input sql function
invbid=drop_na(invbid)

#search by Pedigree
peds=data.frame(id=unique(c(dhSP$Female_pedigree,dhSP$`Male tester_pedigree`)))# input sql function
peds=drop_na(peds)

#### run sql ####
#check AWS results! run once
source('/repos/ELBFA_GS_modules/geno_mod/func_sql.R')
#sql_origin(miss_inv3=invbid,season=season,year=year,method="barcode")
#sql_origin(miss_inv3=peds,season=season,year=year,method="peds")

s3load(object = 'ycao1/ITG/Cucumber/sql_origin/sql_invbid_name_summer2020.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
s3load(object = 'ycao1/ITG/Cucumber/sql_origin/sql_pedigree_name_summer2020.Rdata', 
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

#### load genotypes
# these genotypes were the result of the inbred_format function from the inbred_geno_format.R

#load("/mnt/DH_processing/DH_2020_geno_cucumber.Rdata")#2020
#s3save(DH_2020, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = 'ycao1/ITG/Cucumber/DH_2020_geno_cucumber.Rdata')

#load("/mnt/inbred_geno_cucumber.Rdata")#training
#s3save(inbred_geno, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = 'ycao1/ITG/Cucumber/inbred_geno_cucumber.Rdata')

#load("/mnt/DH_processing/summer_DH_geno_cucumber.Rdata")#2021
#s3save(summer_DH, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = 'ycao1/ITG/Cucumber/summer_DH_geno_cucumber.Rdata')

s3load(object = 'ycao1/ITG/Cucumber/inbreds/DH_2020_geno_cucumber.Rdata', #2020
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
s3load(object = 'ycao1/ITG/Cucumber/inbreds/inbred_geno_cucumber.Rdata', #training 2015-2019
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')
s3load(object = 'ycao1/ITG/Cucumber/inbreds/summer_DH_geno_cucumber.Rdata', #2021
       bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')


#### Geno format
inbred_geno=dplyr::rename(inbred_geno,c("ProgenyGermID"="GermID") )
seas=get(paste0(season,"_DH"))
inb=rbind(seas,DH_2020,inbred_geno)

################################
#### table pairwise parents ####
################################
#@Generates a table with hybrid pedigree name,male and female invbid and male and female GermID
#the Germ id is needed to match with the genotype file
tst=ids_par(miss_inv5=invbid,inbredAll=Allinbred3,inbredAll2=Allinbred1,inb=inb,dhSP=dhSP1,dhSP2=dhSP2,year="2020",season="summer")

#########################
#### geno processing ####
#########################
#@Generates the synth hybrids using the genotypes in HD format "raw" provided as is by the genomics team, uses as a reference the tst table 
### combine all genotype data use as geno reference
### load genotypes in HD form for season/year to create hybrids
### use any marker data reference for marker name and number
markerdata <- s3read_using(fread, 
                           object = 'Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/MarkerMetadata.txt',
                           bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

#DH_2020_geno_cucumber
geno_DH <- s3read_using(fread, 
                        object = "Cucumber/LongDutch/Y3/2500Samples_PMC_TrainingSet_20211026/HDgenotypes.txt.gz", 
                        bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')


#summer_DH_geno_cucumber.Rdata
geno_DH_1 <- s3read_using(fread, 
                          object = "Cucumber/LongDutch/Y3/DH0_Fall2021/HDgenotypes.txt.gz", 
                          bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')


#inbred_geno_cucumber.Rdata
geno_DH_2 <- s3read_using(fread, 
                          object = "Cucumber/LongDutch/Y3/SeasonSetTrainingSet_20210825_Auto/HDgenotypes.txt.gz", 
                          bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

geno_DH=data.frame(geno_DH)
geno_DH_1=data.frame(geno_DH_1)
geno_DH_2=data.frame(geno_DH_2)

####
geno_DH=geno_DH[!as.numeric(geno_DH$ProgenyGermID)%in%as.numeric(geno_DH_1$ProgenyGermID),]
all_geno_DH =rbind(geno_DH,geno_DH_1)
geno_DH1=all_geno_DH[!as.numeric(all_geno_DH$ProgenyGermID)%in%as.numeric(geno_DH_2$ProgenyGermID),]
all_geno_DH1 =rbind(geno_DH1,geno_DH_2)

#### function geno_proc parameters ####
Dh_data<-tst
hybridsWithBothParentsGT <- tst%>%distinct(pedigree)
finalFileName <- paste0("/mnt/Y3_",season, "SeasonDHSet1",year,"_SynthHybrids_tst.tab")
geno_proc(Dh_data=tst,all_geno_DH1=all_geno_DH1,finalFileName=finalFileName)

#####################
#### dosage file ####
#####################
#@formats the synth hybrids file to the 0,1,2, results can be uploaded to AWS and later combined with inbred data for prediction
file=finalFileName
crop='Cucumber'
PCM=geno_format(file,crop)
#s3save(PCM, bucket = "veg-apd-sdi-predictiveanalytcs-prod-workspace", object = 'ycao1/ITG/Cucumber/DH_2020_geno_cucumber.Rdata')


  