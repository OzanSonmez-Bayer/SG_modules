obj <-get_object("s3://veg-apd-sdi-predictiveanalytics-prod-geno-data/Cucumber/LongDutch/Y3/SeasonSetTrainingSet_20210825_Auto/SeasonSetSynthHybrids/Y3_SeasonSet_SynthHybrids.tab")  
csvcharobj <- rawToChar(obj)  
con <- textConnection(csvcharobj)  
hybrid_geno <- read.csv(file = con, sep="\t", header=T)
unique(hybrid_geno$HybridGermID)

hybrid_wide <- hybrid_geno %>% 
  dplyr::select(HybridGermID, MRN, Genotype) %>% 
  dcast(HybridGermID ~ MRN, value.var="Genotype")

# Turn MRN to 0, 1, 2
hybrid_wide_012 <- apply(hybrid_wide[, 2:ncol(hybrid_wide)], MARGIN = 2, 
                         FUN = function(x){
                           unlist(lapply(x, function(y) sum(as.numeric(unlist(strsplit(y, '|'))[c(1,3)]))))
                         })

hybrid_wide_012 <- cbind(data.frame("HybridGermID"=hybrid_wide[,1]), hybrid_wide_012)

hybrid_wide_ped <- hybrid_wide_012 %>% 
  mutate(HybridGermID = as.numeric(HybridGermID)) %>% 
  dplyr::inner_join(CropIDs_sub, by = c('HybridGermID' = 'ProgenyGermID'))


#### Format inbred and hybrid data and save in a combined file for analyses

inbred_geno <-geno_wide_ped

colnames(inbred_geno)[which(colnames(inbred_geno) == 'ProgenyGermID')] <- 'GermID'
colnames(hybrid_wide_ped)[which(colnames(hybrid_wide_ped) == 'HybridGermID')] <- 'GermID'


hybrid_wide_ped <- hybrid_wide_ped %>% 
  select(GermID, Pedigree, colnames(hybrid_wide_ped)[which(!colnames(hybrid_wide_ped) %in% c('GermID', 'Pedigree'))]) %>% 
  mutate(GermID = as.numeric(GermID))

# Check to see if there are any NAs in hybrid geno file

tst <- hybrid_wide_ped[, 3:ncol(hybrid_wide_ped)]
tst <- as.matrix(tst)
for (i in 1:ncol(tst)){
  tst[,i][which(is.na(tst[,i]))] <- ceiling(mean(tst[,i], na.rm = T))
}

tst <- as.data.frame(tst)
hybrid_wide_ped[,3:ncol(hybrid_wide_ped)] <- tst
dim(hybrid_wide_ped)

save(hybrid_wide_ped,file="hybrid_geno_cucumber.Rdata")
