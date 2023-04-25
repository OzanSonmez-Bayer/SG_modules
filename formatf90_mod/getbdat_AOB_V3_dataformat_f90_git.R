# get cleaned bdat based on dat_2
getbdat_f90 = function(dataset, n_gen=5){
  #dataworker 
  #bmrd.field_name:bmrd.rep_number 
  #bmrd.field_name 
  #bmrd.rep_number 
  #bmrd.grow_year 
  #bmrd.yr_loc 
  #ped(bmrd.pedigree_name) 
  #R 
  
  #ITG 
  #BR_FIELD_ID:TREP 
  #BR_FIELD_ID 
  #TREP 
  #GROWSEASON 
  #REPETITION 
  #ped(PEDIGREE_NAME) 
  #R 
  
  if (is.null(CropIDs)){ # send warning if crop_ids is NULL, i.e. didn't load
    "A-matrix BLUP model cannot run without loading s3 bucket CropIDs dataframe. Please check and try again."
  }
  
  #dataset<-dat
  #trtList <- unique(dataset$OBSRVTN_REF_CD)
  #trait   <- trtList[1]
  
  bdat <- dataset %>% # convert to correct types across dataframe for trait data
    mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
    mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
    #mutate(bmrd.grow_month = as.character(bmrd.grow_month)) %>% 
    mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME)) %>% 
    mutate(FIELD_NAME = as.character(FIELD_NAME)) %>% 
    mutate(ORIGIN = as.character(ORIGIN)) %>%
    mutate(GERMPLASM_ID = as.character(as.numeric(GERMPLASM_ID)))%>%
           mutate(pedigree = PEDIGREE_NAME)
  
  print("starting making ped file")
  
  peds = make_ped_file(df = bdat, crop_ids=CropIDs, n_gen = n_gen)
  
  print("finish making ped file")
  
  # add other components to the model formula
  
  bdat$TREP=paste(bdat$REP_NUMBER,bdat$TRACK_ID,sep='_')
  bdat$TREP<-as.factor(bdat$TREP)
  
  # add other components to the model formula
  if (nrow(bdat) > 100){ # only model a trait if there's at least 100 data points
    cols_chk <- c('FIELD_NAME', 'GROWSEASON', 'REPETITION', 'TREP')
    levs <- sapply(bdat[cols_chk],function(x)length(unique(x)))
    print(levs)
    rando <- ""
    #### replace with bmrd.br_field_id ###############
    if (levs['FIELD_NAME'] > 1){
      rp=sapply(unlist(lapply(split(bdat,bdat$FIELD_NAME,drop=T),
                              function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
      if(sum(rp, na.rm = T) >= 1){
        rando <- '~FIELD_NAME+TREP'
      }else{
        rando <- '~FIELD_NAME'
      }
    }else{
      if (levs['TREP'] > 1){
        rando <- '~TREP'
      }
    }
    if (levs['GROWSEASON'] > 1){
      if (rando != ""){
        rando <- paste0(rando, '+ GROWSEASON') 
      }else{
        rando <- '~GROWSEASON'
      }
    }
    if (levs['REPETITION'] > 1){
      if (rando != ''){
        rando <- paste0(rando, '+ REPETITION')
      }else{
        rando <- '~REPETITION'
      }
    }
    
    if (rando != ""){
      randof <- formula(paste0(rando, '+ped(ID)')) #ped(bmrd.pedigree_name) -> ped(ID)
    }else{
      randof <- formula('~ped(ID)')
    }
    print(randof)
    fix <- formula('TRAIT_VALUE~1')
  }
    
  pd <- peds$ped
  pd = unique(pd)
  
  
  # replace n_gen with Inf inbreedings and calculate A matrix inverse here:
  if (sum(is.infinite(pd$inbreeding))>0){
    pd0 = pd[is.infinite(pd$inbreeding), ]
    pd0dir = paste0('pd_Inf.RData')
    save(pd0, file=pd0dir)
    pd$inbreeding[is.infinite(pd$inbreedin)] = n_gen
  }
  
  

#### library MCMCglmm
  #library(MCMCglmm)
  #pd1<-pd[,1:3]
  #colnames(pd1)<-c("id","dam","sire")
  #pd2<-inverseA(pedigree=pd1, nodes="ALL", scale=TRUE, reduced=FALSE,
  #         tol = .Machine$double.eps^0.5)
  #pd$inbreeding <- pd2$inbreeding
  
  
  #ord <- orderPed(pd)
  #row_no = pd[which(ord==-1),]
  #pd <- pd[order(ord),]
  #pd = pd%>%filter(!ID==row_no$ID)
  
  pd.ainv <- asreml.Ainverse(pd, fgen = c('inbreeding', n_gen))
  inbreeding = pd.ainv$inbreeding
  pd$inbreeding <- inbreeding
  #pd <- data.frame(peds$pd, inbreeding)
  
  pd <- pd[!duplicated(pd[,1]),]
  ##names(pd) <- c('Pedigree', 'P1', 'P2', 'inbreeding')
  
  my_list = list("bdat" = bdat, "random" = randof, "fix" = fix, "pd"=pd)
  return(my_list) 
}

## UOM CONVERSION
conv_uom = function(pheno){
  trait = unique(pheno$OBSRVTN_REF_CD)
  
  pheno$UOM=as.character(pheno$UOM) 
  kv = data.frame(fts=as.character(c('Miles',"Inches","Kilometers","Meters",'Centimeters',"Millimeters","Kilograms","Grams","Pounds","Ounces","Grams ENG","Centimeters ENG","Kilograms/Plot","Pounds/Plot","Ounces/Plot")),
                  conv = as.character(c("mi","inch","km",'m',"cm","mm",'kg','g',"lbs","oz","g","cm","kg","lbs","oz")))
  for (i in 1:length(trait)){
    
    dat=subset(pheno,OBSRVTN_REF_CD==trait[i])
    tab=sort(table(dat$UOM))
    tab=tab[tab>0]
    if (length(tab) > 1){
      tmp1=mean(as.numeric(subset(dat,UOM==as.character(names(tab)[1]))$TRAIT_VALUE))
      tmp2=mean(as.numeric(subset(dat,UOM==as.character(names(tab)[2]))$TRAIT_VALUE))
      if (!is.na(tmp1) & !is.na(tmp2)){
        if (tmp1 > tmp2){
          rto = round(tmp1/tmp2)
        }else{
          rto = round(tmp2/tmp1)
        }
        if (rto >= 2){
          if (as.character(names(tab)[1]) %in% kv$fts & as.character(names(tab)[2]) %in% kv$fts){
            uom1=as.character(subset(kv,fts==as.character(names(tab)[1]))$conv)
            uom2=as.character(subset(kv,fts==as.character(names(tab)[2]))$conv)
            conv=round(conv_unit(as.numeric(subset(dat,UOM==as.character(names(tab)[1]))$TRAIT_VALUE),uom1,uom2),2)
            uom_cnv <- unique(dat$TRAIT_VALUE_uom_id[dat$UOM == as.character(names(tab)[2])])
          }else{
            #Generic function
            dat$TRAIT_VALUE=as.numeric(dat$TRAIT_VALUE)
            mn1=mean(as.numeric(dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE),na.rm=T)
            mn2=mean(as.numeric(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE),na.rm=T)
            rtio=mn1/mn2
            if (round(rtio)>2) {
              conv=round(dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE/rtio,3)
              uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE_uom_id)
            }else{
              conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
              uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE_uom_id)
            }
          } # END IF UNITS NOT FOUND
          
        }else{
          conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
          uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE_uom_id)
        } # END IF RATION OF UNITS > 3
      }else{
        conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
        uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE_uom_id)
      } # end if tmp 1 or tmp 2 is NA
    }else{
      conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
      uom_cnv <-  unique(dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE_uom_id)
    } # END IF THERE'S ONLY ONE UNIT
    pheno[pheno$OBSRVTN_REF_CD==trait[i] & pheno$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE=conv
    pheno[pheno$OBSRVTN_REF_CD==trait[i] & pheno$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE_uom_id <- uom_cnv
  }
  pheno
}




