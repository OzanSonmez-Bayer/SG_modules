
getoutput_f90 = function(out,dataset,trait,gg,n_gen=5){
  #@out : folder name of the analysis
  #@dataset : phenotypic file after UOM conversion
  #@trait : input trait from traitlist AOB onboarding
  #@gg : genotypic imputed markers file, if NULL will run pedigree P-BLUP
  #@n_gen: number of generations back used in deep pedigree 


  outbox = paste0(out, "SSGBLUP",'_deep/', trait, '_deep')
  if(file.exists(outbox)==F) {dir.create(outbox, recursive = T)}
  
  setwd(outbox)
  dat1  = dataset[is.element(dataset$OBSRVTN_REF_CD, trait), ]
  
  ####################
  
  # add P1 and P2
  my_list = getbdat_f90(dataset=dat1, n_gen=5)
  
  bdat    = my_list$bdat
  random0 = my_list$random # save for future blupf90 models
  fix0    = my_list$fix
  pd      = my_list$pd
  
  print("bdat,random effects and pedigree finished")
  ##########
  # pedigree initial formatting
  ##########
  
  tmp = left_join(pd, cropids, by=c('ID'='M.GERMPLASM.X_ID'))
  tmp = left_join(tmp, cropids, by=c('PARENT_FEMALE'='M.GERMPLASM.X_ID'))
  tmp = left_join(tmp, cropids, by=c('PARENT_MALE'='M.GERMPLASM.X_ID'))
  tmp = tmp[, c('M.GERMPLASM.PEDIGREE.x', 'M.GERMPLASM.PEDIGREE.y', 'M.GERMPLASM.PEDIGREE')]
  colnames(tmp) =  c("ID", "PARENT_FEMALE", "PARENT_MALE")
  tmp$inbreeding =pd$inbreeding
  names(tmp) <- c('Pedigree', 'P1', 'P2', 'inbreeding')
  pd<-tmp%>%dplyr::select('Pedigree', 'P1', 'P2', 'inbreeding')
  
  # simplify pedigree ID 
  
  if(is.null(gg) == FALSE){
    id = data.frame('PEDIGREE_NAME'=c(as.character(bdat$PEDIGREE_NAME), as.character(pd$Pedigree), as.character(pd$P1), as.character(pd$P2), as.character(gg$PEDIGREE)), 'PEDIGREE_NAME_ID'=NA)
  }else{
    id = data.frame('PEDIGREE_NAME'=c(as.character(bdat$PEDIGREE_NAME), as.character(pd$Pedigree), as.character(pd$P1), as.character(pd$P2)), 'PEDIGREE_NAME_ID'=NA)
  }
  
  
  id = id[duplicated(id)==F, ]
  id$PEDIGREE_NAME    = as.factor(id$PEDIGREE_NAME)
  id$PEDIGREE_NAME_ID = as.numeric(id$PEDIGREE_NAME)
  id$PEDIGREE_NAME_ID = str_pad(id$PEDIGREE_NAME_ID, width=max(nchar(id$PEDIGREE_NAME_ID),na.rm=T), side="left", pad="0")
  id$PEDIGREE_NAME_ID = paste0('ID_', id$PEDIGREE_NAME_ID)
  id$PEDIGREE_NAME_ID[is.na(id$PEDIGREE_NAME)] = 0
  
  #pd0 = pd
  pd  = left_join(pd, id, by=c('Pedigree'='PEDIGREE_NAME'))
  pd  = left_join(pd, id, by=c('P1'='PEDIGREE_NAME'))
  pd  = left_join(pd, id, by=c('P2'='PEDIGREE_NAME'))
  pd  = pd[, c('PEDIGREE_NAME_ID.x', 'PEDIGREE_NAME_ID', 'PEDIGREE_NAME_ID.y')]
  colnames(pd) = c('PEDIGREE_NAME', 'P2', 'P1')
  pd$V4 = 0
  pd$V5 = 0
  pddir = paste0(outbox, '/rawpedigree.txt')
  write.table(pd, file=pddir, quote=F, sep=' ', row.names=F, col.names =F)
  
  #### equivalency ID and pedigree name needed to output GEBV 
  #### for every ID code and not only for ID with NA in phenotype file
  outids<-data.frame(id)
  outidsdir = paste0(outbox, '/outids.txt')
  #write.table(outids, file=outidsdir, quote=F, sep=',', row.names=F, col.names =F)
  
  
  ##########
  # genotype initial formatting
  ##########
  if(is.null(gg)==FALSE){
  s = gg
  s = s[is.element(s$PEDIGREE, c(bdat$PEDIGREE_NAME, bdat$P1, bdat$P2)), ]
  s = left_join(s, id, by=c('PEDIGREE'='PEDIGREE_NAME'))
  PEDIGREE_NAME_ID = s$PEDIGREE_NAME_ID
  s<-data.frame(s)
  ss = s[, is.element(colnames(s), c('PEDIGREE', 'ProgenyGermID', 'PEDIGREE_NAME_ID'))==F]
  ss = as.data.frame(ss)
  colnames(ss)<-paste0("m",seq(1:ncol(ss)))
  dat<-data.frame("ID"=PEDIGREE_NAME_ID,ss)
  }
  #########
  # random effects
  #########
  
  pool = unlist(strsplit(as.character(random0), "\\~|\\+| "))
  pool = pool[is.element(pool, "")==F] 
  random = intersect(pool, colnames(bdat))
  if (fix0=='TRAIT_VALUE ~ 1') {bdat$intercept = 1}
  
  ###########
  # format input data
  ###########
  
  #### rawdata JWAS
  
  raw = bdat[, c('PEDIGREE_NAME','TRAIT_VALUE', 'intercept', random)]
  raw = left_join(raw, id)
  str(raw$PEDIGREE_NAME)
  rawdir = paste0(outbox, '/rawdata.txt')
  
  #write.table(raw, file= rawdir, quote=F, sep=' ', row.names=F, col.names =F)
  
  ### 
  ### exception if crossvalidation is requested
  
  #if (is.null(seas)== TRUE){
  #  raw=raw
  #}else{
  #  test <- raw %>% filter(GROWSEASON == seas) %>% .$PEDIGREE_NAME %>% unique()
  #  raw <- raw %>%
  #    mutate(TRAIT_VALUE = ifelse(PEDIGREE_NAME %in% test, NA, TRAIT_VALUE)) 
  #}
  
  #### phenos JWAS
  phenos_raw<-raw[,colnames(raw)%in%c('PEDIGREE_NAME_ID','TRAIT_VALUE', 'intercept')]
  colnames(phenos_raw)<-c("y","intercept","ID")
  if(c("TREP")%in%colnames(raw)==TRUE){levels(raw$TREP)<-LETTERS[seq(1:nlevels(raw$TREP))]}
  if(c("GROWSEASON")%in%colnames(raw)==TRUE){raw$GROWSEASON<-as.factor(raw$GROWSEASON)}
  if(c("GROWSEASON")%in%colnames(raw)==TRUE){levels(raw$GROWSEASON)<-LETTERS[seq(1:nlevels(raw$GROWSEASON))]}
  if(c("REPETITION")%in%colnames(raw)==TRUE){raw$REPETITION<-as.factor(raw$REPETITION)}
  if(c("REPETITION")%in%colnames(raw)==TRUE){levels(raw$REPETITION)<-LETTERS[seq(1:nlevels(raw$REPETITION))]}
  phenos_raw2<-raw[,!colnames(raw)%in%c('PEDIGREE_NAME_ID','TRAIT_VALUE', 'intercept','PEDIGREE_NAME')]
  colnames(phenos_raw2)<-paste0("x",1:ncol(phenos_raw2))
  phenos_jwas<-data.frame("ID"= phenos_raw$ID,phenos_raw2,"y"= phenos_raw$y,"intercept"=phenos_raw$intercept)
  jlphenodir = paste0(outbox, '/testphenos_jwas.txt')
  write.table(phenos_jwas, file=jlphenodir, quote=F, sep=',', row.names=F, col.names =T)
  
  print("writing phenotype file")
  
  #### pedigree JWAS
  
  pedigree<-read.table(paste0(outbox,'/rawpedigree.txt'))
  ped_jwas<-data.frame("ID"=pedigree$V1,"Sire"=pedigree$V2,"Dam"=pedigree$V3)
  parents1<-data.frame("ID"=ped_jwas$Sire,"Sire"=0,"Dam"=0) 
  parents2<-data.frame("ID"=ped_jwas$Dam,"Sire"=0,"Dam"=0) 
  all_ped<-rbind(ped_jwas,parents1,parents2)
  all_ped=all_ped[!duplicated(all_ped$ID),]
  jlpeddir = paste0(outbox, '/testped_jwas.txt')
  write.table(all_ped, file=jlpeddir, quote=F, sep=',', row.names=F, col.names =T)
  
  outebv<-data.frame(all_ped$ID)
  outebvdir = paste0(outbox, '/outebv.txt')
  #write.table(outebv, file=outebvdir, quote=F, sep=',', row.names=F, col.names =F)
  
  print("writing pedigree file")
  
  #### genotypes JWAS
  #library(fread)
  
  if(is.null(gg)==FALSE){
  genos<-dat[dat$ID%in%all_ped$ID,]
  jlgenodir = paste0(outbox, '/testgenos_jwas.txt')
  #fwrite
  write.table(genos,file=jlgenodir,quote=F, sep=',', row.names=F, col.names =TRUE)
  print("writing genotype file")
  }
  
 
}
