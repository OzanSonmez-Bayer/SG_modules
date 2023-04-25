# as.fun
as.fun <- function(x,rando,mdl,fix,pd.inv=NA,fm=asreml.gaussian()){
  if(mdl == 'blu'){
    fix=formula('TRAIT_VALUE~CROSS_NAME')
    if (rando == ""){
      as_mdl <- do.call(asreml, args=list(fixed=fix,data=x,
                                          na.method.X='include',trace=FALSE,family=fm,
                                          workspace = 320e+6, pworkspace = 320e+6)
      )
    }else{
      as_mdl <- do.call(asreml, args=list(fixed=fix,random=formula(rando),data=x,
                                          na.method.X='include',trace=FALSE,family=fm,
                                          workspace = 320e+6, pworkspace = 320e+6)
      )
    }
    # as_mdl <- do.call(asreml, args=list(fixed=fix,random=rando,data=x,
    #                                     na.method.X='include',trace=FALSE,family=fm)
    # )
  }
  return(as_mdl)
}



# BLUP functions
pedSingle <- function(ped){
  if (length(gregexpr(')',ped)[[1]]) >1){
    x=unlist(strsplit(ped,')',fixed=T))
    y=gregexpr("\\.[A-z0-9][A-z0-9]|>[A-z0-9][A-z0-9]", x[length(x)])[[1]]
    if (length(y)>3){
      cped=paste(paste(x[-length(x)],collapse=')'),substr(x[length(x)],1,y[4]),collapse ='')
    }else{cped=ped}
  }else{
    split.str= grep("\\.[A-z0-9][A-z0-9]|>[A-z0-9][A-z0-9]", ped)
    split.pos <- gregexpr("\\.[A-z0-9][A-z0-9]|>[A-z0-9][A-z0-9]", ped)[[1]]
    if(length(split.pos)>3){
      split.length <- attr(split.pos, "match.length")
      split.start <- sort(c(split.pos, split.pos+split.length))
      split.end <- c(split.start[-1]-1, nchar(ped))
      y=paste(substring(ped,split.start,split.end)[1:6],collapse='')     
      cped=paste(substr(ped,1,split.pos[1]-1),y,collapse='')
    }else{cped=ped}
  }
  cped
}


RunBLUE <- function(bdat, listofcrosses) {
  
  bdat$TREP <- paste(bdat$REP_NUMBER,bdat$TRACK_ID,sep='_')
  
  print("Running blue")
  
  bdat_subset <- subset(bdat, CROSS_NAME %in% listofcrosses)
  
  ###Recalculate rando
  cols=c('GROWSEASON','BR_LOC_ID','PLOT_ID','ABS_R','ABS_C','TEST_SET_NAME','SET_OWNER_PROG','BR_FIELD_ID','FIELD_NAME','BR_REP_ID','REP_NUMBER','TREP','COUNTRY','SUB_COUNTRY', 'REPETITION')
  levs=sapply(bdat[cols],function(x)length(levels(x)))
  print(levs)
  rando=''
  if (levs['GROWSEASON'] >1) {
    #Looking for logic to account for Growseason. Appears to be reaching boundary 
    #ss=sapply(unlist(lapply(split(bdat,bdat$GROWSEASON,drop=T),function(x)unique(x$CROSS_NAME))),function(x) sum(x) > 1))
    rando='~GROWSEASON'
  }
  if (!is.na(levs['REPETITION'])){
    if (levs['REPETITION'] > 1){
      if (rando != ''){
        rando <- paste0(rando, '+REPETITION')
      }else{
        rando <- '~REPETITION'
      }
    }
  }
  
  # if (levs['FIELD_NAME'] >1) {
  #   #Add rep_number as covariate? Only if reps within same field.
  #   rp=sapply(unlist(lapply(split(bdat,bdat$FIELD_NAME,drop=T),function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
  #   if (rando != ''){
  #     #If multiple reps per field
  #     if (sum(rp,na.rm=T) >= 1){
  #       rando=paste0(rando,'+FIELD_NAME:TREP')
  #     }else{
  #       rando=paste0(rando,'+FIELD_NAME')
  #     }
  #   }else{
  #     if (sum(rp,na.rm=T) >= 1){
  #       rando=paste0(rando,'~FIELD_NAME:TREP')
  #     }else{
  #       rando='~FIELD_NAME'
  #     }
  #   }
  # }else{
  #   #if (levs['REP_NUMBER']>1){
  #   if (levs['TREP']>1){
  #     if (rando != ''){
  #       #rando=paste0(rando,'+REP_NUMBER')
  #       rando=paste0(rando,'+TREP')
  #     }else{
  #       rando=('~TREP') 
  #     }
  #   }
  #   # rm(rp)
  # }
  
  if (levs['BR_FIELD_ID'] >1) {
    #Add rep_number as covariate? Only if reps within same field.
    rp=sapply(unlist(lapply(split(bdat,bdat$FIELD_NAME,drop=T),function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
    if (rando != ''){
      #If multiple reps per field
      if (sum(rp,na.rm=T) >= 1){
        rando=paste0(rando,'+BR_FIELD_ID:TREP')
      }else{
        rando=paste0(rando,'+BR_FIELD_ID')
      }
    }else{
      if (sum(rp,na.rm=T) >= 1){
        rando=paste0(rando,'~BR_FIELD_ID:TREP')
      }else{
        rando='~BR_FIELD_ID'
      }
    }
  }else{
    #if (levs['REP_NUMBER']>1){
    if (levs['TREP']>1){
      if (rando != ''){
        #rando=paste0(rando,'+REP_NUMBER')
        rando=paste0(rando,'+TREP')
      }else{
        rando=('~TREP') 
      }
    }
    # rm(rp)
  }
  
  if (bdat$OBSRVTN_REF_CD[1]  == 'YIELDWT_HAARE'){
    if (rando != ''){
      rando <- paste0(rando, '+STAND_HAARE')
    }else{
      rando <- '~STAND_HAARE'
    }
  }
  rm(levs)
  ###
  
  pred_blu=data.frame()
  as_blu <- as.fun(bdat_subset, rando, "blu")
  
  
  pred_blu <- as.data.frame(predict(as_blu,classify = 'CROSS_NAME', max.iter=1)$pred$pval)
  
  pred_blu
  
}


blue_fun <- function(dataset_1, crop, listofcross){
  
  dataset <- subset(dataset_1, CROSS_NAME %in% listofcrosses)
  
  outDF <- data.frame()
  trtList <- unique(dataset$OBSRVTN_REF_CD)
  dataset$TREP <- paste(dataset$REP_NUMBER,dataset$TRACK_ID,sep='_')
  for (i in 1:length(trtList)){
    # Some traits are not in the dataset
    # Check the existence of trait
    trait <- trtList[i]
    print(trait)
    bdat <- dataset %>% 
      filter(OBSRVTN_REF_CD == trait) %>% 
      mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
      mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
      mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME)) %>% 
      mutate(BR_FIELD_ID = as.character(BR_FIELD_ID)) %>% 
      mutate(ORIGIN = as.character(ORIGIN))
    
    
    if(sum(bdat$ORIGIN == '')>0){
      bdat$ORIGIN[which(bdat$ORIGIN == '')] <- bdat$PEDIGREE_NAME[which(bdat$ORIGIN == '')]
    }
    
    if (!'P1' %in% colnames(bdat)){
      print("parentage pedigree doesn't exit")
      bdat <- bdat %>% 
        separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
      p1 <- unlist(lapply(bdat$P1, FUN = pedSingle))
      p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      bdat$P1 <- p1
      bdat$P2[bdat$P2=='']=NA
      p2 <- unlist(lapply(bdat$P2, FUN = pedSingle))
      p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      bdat$P2 <- p2
      
    }else{
      print("parentage pedigree exist")
      colnames(bdat)[which(colnames(bdat) == 'parent.female_pedigree_name')] = 'P1'
      colnames(bdat)[which(colnames(bdat) == 'parent.male_pedigree_name')] = 'P2'
      
    }
    if (sum(is.na(bdat$P1)) >0) {
      bdat[is.na(bdat$P1), 'P2'] <- NA
    }
    
    if (sum(is.na(bdat$P2)) >0 ){
      bdat[is.na(bdat$P2), 'P1'] <- NA
    }
    
    bdat$P1 <- as.character(bdat$P1)
    bdat$P2 <- as.character(bdat$P2)
    print(dim(bdat))
    if (nrow(bdat) > 100){
      cols_chk <- c('BR_FIELD_ID', 'GROWSEASON', 'REPETITION', 'TREP')
      bdat[cols_chk] <- lapply(bdat[cols_chk], factor)
      levs <- sapply(bdat[cols_chk],function(x)length(unique(x)))
      print(levs)
      rando <- ""
      if (levs['BR_FIELD_ID'] > 1){
        rp <- sapply(unlist(lapply(split(bdat,bdat$BR_FIELD_ID,drop=T),function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
        if(sum(rp, na.rm = T) >= 1){
          rando <- '~BR_FIELD_ID:TREP'
        }else{
          rando <- '~BR_FIELD_ID'
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
      # 
      # if (levs['TREP'] > 1){
      #   if (rando != ''){
      #     rando <- paste0(rando, '+ TREP')
      #   }else{
      #     rando <- '~ TREP'
      #   }
      # }
      
      
      fix <- formula('TRAIT_VALUE~CROSS_NAME')
      

      
      # pd.aniv <- asreml.Ainverse(pd, fgen = c('Selfing',5))
      # pd.aniv <- asreml.Ainverse(pd)
      
      gca_mdl <- do.call(asreml, args = list(fixed = fix, random = formula(rando), data = bdat, 
                                             na.method.X = 'include', 
                                             trace = FALSE, 
                                             na.method.Y = 'include',
                                             workspace = 320e+6, pworkspace = 320e+6))
      
      pred_df <- as.data.frame(predict(gca_mdl, classify = 'CROSS_NAME', max.iter = 1)$predictions$pvals)
      
      return(pred_df)
    }
  }
}


Ablup_fun <- function(dataset, n_gen, CropIDs, crop){
  outDF <- data.frame()
  trtList <- unique(dataset$OBSRVTN_REF_CD[!is.na(dataset$OBSRVTN_REF_CD)])
  dataset$TREP <- paste(dataset$REP_NUMBER,dataset$TRACK_ID,sep='_')
  for (i in 1:length(trtList)){
    # Some traits are not in the dataset
    # Check the existence of trait
    trait <- trtList[i]
    print(trait)
    bdat <- dataset %>% 
      filter(OBSRVTN_REF_CD == trait) %>% 
      mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
      mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
      mutate(pedigree = as.character(PEDIGREE_NAME)) %>% 
      mutate(BR_FIELD_ID = as.character(BR_FIELD_ID)) %>% 
      mutate(ORIGIN = as.character(ORIGIN))
    
    if(sum(bdat$ORIGIN == '')>0){
      bdat$ORIGIN[which(bdat$ORIGIN == '')] <- bdat$pedigree[which(bdat$ORIGIN == '')]
    }
    
    if (!'P1' %in% colnames(bdat)){
      print("parentage pedigree doesn't exit")
      bdat <- bdat %>% 
        separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
      p1 <- unlist(lapply(bdat$P1, FUN = pedSingle))
      p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      bdat$P1 <- p1
      bdat$P2[bdat$P2=='']=NA
      p2 <- unlist(lapply(bdat$P2, FUN = pedSingle))
      p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      bdat$P2 <- p2
      
    }else{
      print("parentage pedigree exist")
      colnames(bdat)[which(colnames(bdat) == 'parent.female_pedigree_name')] = 'P1'
      colnames(bdat)[which(colnames(bdat) == 'parent.male_pedigree_name')] = 'P2'
      
    }
    if (sum(is.na(bdat$P1)) >0) {
      bdat[is.na(bdat$P1), 'P2'] <- NA
    }
    
    if (sum(is.na(bdat$P2)) >0 ){
      bdat[is.na(bdat$P2), 'P1'] <- NA
    }
    
    bdat$P1 <- as.character(bdat$P1)
    bdat$P2 <- as.character(bdat$P2)
    
    if (nrow(bdat) > 50){
      cols_chk <- c('BR_FIELD_ID', 'GROWSEASON', 'REPETITION', 'TREP')
      bdat[cols_chk] <- lapply(bdat[cols_chk], factor)
      levs <- sapply(bdat[cols_chk],function(x)length(unique(x)))
      print(levs)
      rando <- ""
      if (levs['BR_FIELD_ID'] > 1){
        rp <- sapply(unlist(lapply(split(bdat,bdat$BR_FIELD_ID,drop=T),function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
        if(sum(rp, na.rm = T) >= 1){
          rando <- '~BR_FIELD_ID:TREP'
        }else{
          rando <- '~BR_FIELD_ID'
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
      # 
      # if (levs['TREP'] > 1){
      #   if (rando != ''){
      #     rando <- paste0(rando, '+ TREP')
      #   }else{
      #     rando <- '~ TREP'
      #   }
      # }
      
      # new_token(fingerprint = F)
      gen5 <- make_ped_file(bdat, crop_ids = CropIDs, n_gen = n_gen)
      
      outbred <- c('Carrot', 'Onion')
      if (crop %in% outbred){
        A_matrix5 <- asreml.Ainverse(gen5$ped)
      } else{
        A_matrix5 <- asreml.Ainverse(gen5$ped, fgen = list("inbreeding", n_gen))
      }
      
      # bdat <- bdat %>% 
      #   dplyr::select(-ID)
      bdat <- left_join(bdat, gen5$distinct_lines[ , c("ID", "germplasm_id", "pedigree")], by = "pedigree")
      print(colnames(bdat))
      
      if (rando != ""){
        randof <- formula(paste0(rando, '+ped(ID)'))
      }else{
        randof <- formula('~ped(ID)')
      }
      print(randof)
      fix <- formula('TRAIT_VALUE~1')
      
      # bdat$ID <- as.factor(bdat$ID)
      gca_mdl <- do.call(asreml, args = list(fixed = fix, random = randof, data = bdat, 
                                             na.method.X = 'include', trace = FALSE, 
                                             ginverse = list(ID = A_matrix5$ginv),
                                             workspace = 320e+6, pworkspace = 320e+6))
      if (gca_mdl$last.message == 'LogLikelihood not converged'){
        gca_mdl <- update.asreml(gca_mdl)
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(ID)', max.iter = 1)$predictions$pvals)
      }else{
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(ID)', max.iter = 1)$predictions$pvals)
      }
      
      vce <- summary(gca_mdl)$varcomp
      ped_pos <- which(rownames(vce) == "ped(ID)!ped")
      if (median(pred_df$predicted.value)==0 | (sum(c('Boundary','?','Singular') %in% vce[-ped_pos,5]) >= 1 & !vce[ped_pos,5] %in% c('Boundary','?','Singular'))) {
        gca_mdl <- rmboundary.asrtests(asrtests(gca_mdl))[[1]]
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(ID)', max.iter = 1)$predictions$pvals)
        vce <- summary(gca_mdl)$varcomp
      }
      
      h2 <- round(vce['ped(ID)!ped',2]/sum(vce[,2]),4)
      ll <- summary(gca_mdl)$loglik
      print(c("h2: ", h2, "ll: ", ll))
      
      int_val <- as.integer(bdat$TRAIT_VALUE)
      cat_check <- sum(bdat$TRAIT_VALUE - int_val,na.rm=T)
      if (max(bdat$TRAIT_VALUE, na.rm = T) > 9){
        cat_check <- 1
      }
      
      if ((sd(pred_df$predicted.value, na.rm = T)==0)|(cat_check == 0 & (min(pred_df$predicted.value, na.rm = T)<1 | (max(bdat$TRAIT_VALUE, na.rm = T) - max(pred_df$predicted.value, na.rm = T)) > 3.5))){
        print('Poisson Distribution')
        ####CHANGING TO POISSON
        print(c('running poisson, cat_check:',cat_check,'min_pred',min(pred_df$predicted.value, na.rm = T),'max_pred:',max(pred_df$predicted.value,na.rm = T),'sd:',sd(pred_df$predicted.value, na.rm = T)))  
        gca_mdl_pois <- do.call(asreml, args = list(fixed = fix, random = randof, data = bdat, 
                                               na.method.X = 'include', 
                                               trace = FALSE,
                                               family = asreml.poisson(),
                                               ginverse = list(ID = A_matrix5$ginv),
                                               workspace = 320e+6, pworkspace = 320e+6))
        print(summary(gca_mdl_pois)$varcomp)
        
        pred_df_pois <- predict(gca_mdl_pois,classify = 'ped(ID)',max.iter=1)$predictions$pvals
        
        if ( (sd(pred_df_pois$predicted.value, na.rm = T)>0.001) & min(pred_df_pois$predicted.value, na.rm = T) >= 1 ) {
          pred_df[,1:3] <- pred_df_pois[,1:3]
          gca_mdl <- gca_mdl_pois
        } ## IF POISS
      }
      
      rnd_ebv <- data.frame(coef(gca_mdl)$random, gca_mdl$vcoeff$random)
      names(rnd_ebv) <- c('BLUP', 'BLUP.se')
      rnd_ebv <- rnd_ebv[grep('ped', rownames(rnd_ebv)),]
      rnd_ebv$pedigree <- gsub('ped(ID)_', '', rownames(rnd_ebv), fixed = TRUE)
      rnd_ebv$BLUP <- round(rnd_ebv$BLUP,4)
      PEV <- rnd_ebv[,2]*gca_mdl$sigma2
      # reliability <- round(1-PEV/vce['ped(PEDIGREE_NAME)!ped',2],4)
      reliability=round(1-PEV/((1 + A_matrix5$inbreeding)*vce["ped(ID)!ped",2]),4)
      accuracy <- sqrt(reliability)
      rnd_ebv$reliability <- reliability
      
      pred_ebv <- merge(pred_df, rnd_ebv, by.x = 'ID', by.y = 'pedigree', all.x = TRUE, all.y = TRUE)
      pred_ebv$trait <- trait
      pred_ebv$n <- NA
      pred_ebv$N <- nrow(bdat)
      pred_ebv$h2 <- h2
      pred_ebv$likelihood <- ll
      pred_ebv$mdl <- 'ped'
      n_obs <- summaryBy(TRAIT_VALUE ~ PEDIGREE_NAME, data = bdat, FUN = function(x) {length(x)})
      colnames(n_obs)[2] <- 'n'
      
      # Add pedigree back to P_SCA and P_GCA tables
      pred_ebv$ID <- as.character(pred_ebv$ID)
      pred_ebv <- CropIDs %>% 
        dplyr::select(M.GERMPLASM.X_ID, M.GERMPLASM.PEDIGREE) %>% 
        unique() %>% 
        mutate(ID = M.GERMPLASM.X_ID) %>% 
        right_join(pred_ebv, by = 'ID') %>% 
        dplyr::select(-c('ID', 'M.GERMPLASM.X_ID'))
      colnames(pred_ebv)[which(colnames(pred_ebv) == 'M.GERMPLASM.PEDIGREE')] <- 'PEDIGREE_NAME'
      
      pred_ebv <- pred_ebv %>% 
        left_join(n_obs, by = 'PEDIGREE_NAME') %>% 
        mutate(n.x = n.y) %>% 
        dplyr::select(-n.y) %>% 
        dplyr::rename(n = n.x) %>% 
        mutate(mdl = if_else(as.character(PEDIGREE_NAME) %in% unique(bdat$PEDIGREE_NAME[bdat$LTYPE == 'Hybrid']), 'P_SCA', 'P_GCA'))
      
      freq_parent <- table(c(bdat$P1, bdat$P2))
      unique_p <- unique(c(bdat$P1, bdat$P2))
      unique_p <- unique_p[!is.na(unique_p)]
      parent_hybrid <- data.frame('Parent' = unique_p, 'Number_Crosses' = 0)
      for (u_p in 1:length(unique_p)){
        hybrid_temp <- bdat$PEDIGREE_NAME[which(bdat$P1 == unique_p[u_p] | bdat$P2 == unique_p[u_p])]
        parent_hybrid[u_p, 'Number_Crosses'] <- length(unique(hybrid_temp))
      }
      rownames(parent_hybrid) <- parent_hybrid$Parent
      pred_ebv[pred_ebv$mdl == 'P_GCA', ]$n <- parent_hybrid[as.character(pred_ebv[pred_ebv$mdl == 'P_GCA',]$PEDIGREE_NAME), 'Number_Crosses']
      if (i == 1){
        outDF <- pred_ebv
      }else{
        outDF <- rbind.data.frame(outDF, pred_ebv)
      } 
    }
  }
  GCA_result <- subset(outDF, mdl == 'P_GCA')
  SCA_result <- subset(outDF, mdl == 'P_SCA')
  blup_results <- list(GCA_result, SCA_result)
}


pblup_fun_1 <- function(dataset, ainvse_mat){
  outDF <- data.frame()
  trtList <- unique(dataset$OBSRVTN_REF_CD)
  for (i in 1:length(trtList)){
    # Some traits are not in the dataset
    # Check the existence of trait
    trait <- trtList[i]
    print(trait)
    bdat <- dataset %>% 
      filter(OBSRVTN_REF_CD == trait) %>% 
      mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
      mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
      mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME)) %>% 
      mutate(BR_FIELD_ID = as.character(BR_FIELD_ID)) %>% 
      mutate(ORIGIN = as.character(ORIGIN))

    if (nrow(bdat) > 100){
      cols_chk <- c('BR_FIELD_ID', 'GROWSEASON', 'REPETITION', 'TREP')
      bdat[cols_chk] <- lapply(bdat[cols_chk], factor)
      levs <- sapply(bdat[cols_chk],function(x)length(unique(x)))
      print(levs)
      rando <- ""
      if (levs['BR_FIELD_ID'] > 1){
        rp <- sapply(unlist(lapply(split(bdat,bdat$BR_FIELD_ID,drop=T),function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
        if(sum(rp, na.rm = T) >= 1){
          rando <- '~BR_FIELD_ID:TREP'
        }else{
          rando <- '~BR_FIELD_ID'
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
      # 
      # if (levs['TREP'] > 1){
      #   if (rando != ''){
      #     rando <- paste0(rando, '+ TREP')
      #   }else{
      #     rando <- '~ TREP'
      #   }
      # }
      if (rando != ""){
        randof <- formula(paste0(rando, '+ped(PEDIGREE_NAME)'))
      }else{
        randof <- formula('~ped(PEDIGREE_NAME)')
      }
      print(randof)
      fix <- formula('TRAIT_VALUE~1')
      
      gca_mdl <- do.call(asreml, args = list(fixed = fix, random = randof, data = bdat, 
                                             na.method.X = 'include', trace = FALSE, 
                                             ginverse = list(PEDIGREE_NAME = ainvse_mat),
                                             workspace = 320e+6, pworkspace = 320e+6))
      if (gca_mdl$last.message == 'LogLikelihood not converged'){
        gca_mdl <- update.asreml(gca_mdl)
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(PEDIGREE_NAME)', max.iter = 1)$predictions$pvals)
      }else{
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(PEDIGREE_NAME)', max.iter = 1)$predictions$pvals)
      }
      
      vce <- summary(gca_mdl)$varcomp
      print(vce)
      
      pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(PEDIGREE_NAME)', max.iter = 1)$predictions$pvals)
      ped_pos <- which(rownames(vce) == "ped(PEDIGREE_NAME)!ped")
      if (median(pred_df$predicted.value)==0 | (sum(c('Boundary','?','Singular') %in% vce[-ped_pos,5]) >= 1 & !vce[ped_pos,5] %in% c('Boundary','?','Singular'))) {
        gca_mdl <- rmboundary.asrtests(asrtests(gca_mdl))[[1]]
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(PEDIGREE_NAME)', max.iter = 1)$predictions$pvals)
        vce <- summary(gca_mdl)$varcomp
      }
      
      h2 <- round(vce['ped(PEDIGREE_NAME)!ped',2]/sum(vce[,2]),4)
      ll <- summary(gca_mdl)$loglik
      print(c("h2: ", h2, "ll: ", ll))
      
      int_val <- as.integer(bdat$TRAIT_VALUE)
      cat_check <- sum(bdat$TRAIT_VALUE - int_val,na.rm=T)
      if (max(bdat$TRAIT_VALUE, na.rm = T) > 9){
        cat_check <- 1
      }
      
      if ((sd(pred_df$predicted.value, na.rm = T)==0)|(cat_check == 0 & (min(pred_df$predicted.value, na.rm = T)<1 | (max(bdat$TRAIT_VALUE, na.rm = T) - max(pred_df$predicted.value, na.rm = T)) > 3.5))){
        print('Poisson Distribution')
        ####CHANGING TO POISSON
        print(c('running poisson, cat_check:',cat_check,'min_pred',min(pred_df$predicted.value, na.rm = T),'max_pred:',max(pred_df$predicted.value,na.rm = T),'sd:',sd(pred_df$predicted.value, na.rm = T)))  
        gca_mdl_pois <- do.call(asreml, args = list(fixed = fix, random = randof, data = bdat, 
                                                    na.method.X = 'include', 
                                                    trace = FALSE,
                                                    family = asreml.poisson(),
                                                    ginverse = list(PEDIGREE_NAME = pd.aniv$ginv),
                                                    workspace = 320e+6, pworkspace = 320e+6))
        print(summary(gca_mdl)$varcomp)
        
        pred_df_pois <- predict(gca_mdl,classify = 'ped(PEDIGREE_NAME)',max.iter=1)$predictions$pvals
        
        if ( (sd(pred_df_pois$predicted.value, na.rm = T)>0.001) & min(pred_df_pois$predicted.value, na.rm = T) >= 1 ) {
          pred_df[,1:3] <- pred_df_pois[,1:3]
          gca_mdl <- gca_mdl_pois
        } ## IF POISS
      }
      
      rnd_ebv <- data.frame(coef(gca_mdl)$random, gca_mdl$vcoeff$random)
      names(rnd_ebv) <- c('BLUP', 'BLUP.se')
      rnd_ebv <- rnd_ebv[grep('ped', rownames(rnd_ebv)),]
      rnd_ebv$pedigree <- gsub('ped(PEDIGREE_NAME)_', '', rownames(rnd_ebv), fixed = TRUE)
      rnd_ebv$BLUP <- round(rnd_ebv$BLUP,4)
      PEV <- rnd_ebv[,2]*gca_mdl$sigma2
      # reliability <- round(1-PEV/vce['ped(PEDIGREE_NAME)!ped',2],4)
      reliability=round(1-PEV/((1 + pd.aniv$inbreeding)*vce["ped(PEDIGREE_NAME)!ped",2]),4)
      accuracy <- sqrt(reliability)
      rnd_ebv$reliability <- reliability
      
      pred_ebv <- merge(pred_df, rnd_ebv, by.x = 'PEDIGREE_NAME', by.y = 'pedigree', all.x = TRUE, all.y = TRUE)
      pred_ebv$trait <- trait
      pred_ebv$n <- NA
      pred_ebv$N <- nrow(bdat)
      pred_ebv$h2 <- h2
      pred_ebv$likelihood <- ll
      pred_ebv$mdl <- 'ped'
      n_obs <- summaryBy(TRAIT_VALUE ~ PEDIGREE_NAME, data = bdat, FUN = function(x) {length(x)})[,2]
      pred_ebv[rownames(subset(pred_ebv, as.character(PEDIGREE_NAME) %in% unique(bdat$PEDIGREE_NAME))),]$n <- n_obs
      pred_ebv[rownames(subset(pred_ebv, as.character(PEDIGREE_NAME) %in% unique(bdat$PEDIGREE_NAME))),]$mdl <- 'P_SCA'
      pred_ebv[rownames(subset(pred_ebv, !as.character(PEDIGREE_NAME) %in% unique(bdat$PEDIGREE_NAME))),]$mdl <- 'P_GCA'
      
      freq_parent <- table(c(bdat$P1, bdat$P2))
      unique_p <- unique(c(bdat$P1, bdat$P2))
      unique_p <- unique_p[!is.na(unique_p)]
      parent_hybrid <- data.frame('Parent' = unique_p, 'Number_Crosses' = 0)
      for (u_p in 1:length(unique_p)){
        hybrid_temp <- bdat$PEDIGREE_NAME[which(bdat$P1 == unique_p[u_p] | bdat$P2 == unique_p[u_p])]
        parent_hybrid[u_p, 'Number_Crosses'] <- length(unique(hybrid_temp))
      }
      rownames(parent_hybrid) <- parent_hybrid$Parent
      pred_ebv[pred_ebv$mdl == 'P_GCA', ]$n <- parent_hybrid[as.character(pred_ebv[pred_ebv$mdl == 'P_GCA',]$PEDIGREE_NAME), 'Number_Crosses']
      if (i == 1){
        outDF <- pred_ebv
      }else{
        outDF <- rbind.data.frame(outDF, pred_ebv)
      } 
    }
  }
  GCA_result <- subset(outDF, mdl == 'P_GCA')
  SCA_result <- subset(outDF, mdl == 'P_SCA')
  blup_results <- list(GCA_result, SCA_result)
}




pblup_fun <- function(dataset, crop){
  outDF <- data.frame()
  trtList <- unique(dataset$OBSRVTN_REF_CD)
  dataset$TREP <- paste(dataset$REP_NUMBER,dataset$TRACK_ID,sep='_')
  for (i in 1:length(trtList)){
    # Some traits are not in the dataset
    # Check the existence of trait
    trait <- trtList[i]
    print(trait)
    bdat <- dataset %>% 
      filter(OBSRVTN_REF_CD == trait) %>% 
      mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
      mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
      mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME)) %>% 
      mutate(BR_FIELD_ID = as.character(BR_FIELD_ID)) %>% 
      mutate(ORIGIN = as.character(ORIGIN))
    
    
    if(sum(bdat$ORIGIN == '')>0){
      bdat$ORIGIN[which(bdat$ORIGIN == '')] <- bdat$PEDIGREE_NAME[which(bdat$ORIGIN == '')]
    }
    
    if (!'P1' %in% colnames(bdat)){
      print("parentage pedigree doesn't exit")
      bdat <- bdat %>% 
        separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
      p1 <- unlist(lapply(bdat$P1, FUN = pedSingle))
      p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      bdat$P1 <- p1
      bdat$P2[bdat$P2=='']=NA
      p2 <- unlist(lapply(bdat$P2, FUN = pedSingle))
      p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
      bdat$P2 <- p2
      
    }else{
      print("parentage pedigree exist")
      colnames(bdat)[which(colnames(bdat) == 'parent.female_pedigree_name')] = 'P1'
      colnames(bdat)[which(colnames(bdat) == 'parent.male_pedigree_name')] = 'P2'
      
    }
    if (sum(is.na(bdat$P1)) >0) {
      bdat[is.na(bdat$P1), 'P2'] <- NA
    }
    
    if (sum(is.na(bdat$P2)) >0 ){
      bdat[is.na(bdat$P2), 'P1'] <- NA
    }
    
    bdat$P1 <- as.character(bdat$P1)
    bdat$P2 <- as.character(bdat$P2)
    print(dim(bdat))
    if (nrow(bdat) > 100){
      cols_chk <- c('BR_FIELD_ID', 'GROWSEASON', 'REPETITION', 'TREP')
      bdat[cols_chk] <- lapply(bdat[cols_chk], factor)
      levs <- sapply(bdat[cols_chk],function(x)length(unique(x)))
      print(levs)
      rando <- ""
      if (levs['BR_FIELD_ID'] > 1){
        rp <- sapply(unlist(lapply(split(bdat,bdat$BR_FIELD_ID,drop=T),function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
        if(sum(rp, na.rm = T) >= 1){
          rando <- '~BR_FIELD_ID:TREP'
        }else{
          rando <- '~BR_FIELD_ID'
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
      # 
      # if (levs['TREP'] > 1){
      #   if (rando != ''){
      #     rando <- paste0(rando, '+ TREP')
      #   }else{
      #     rando <- '~ TREP'
      #   }
      # }
      if (rando != ""){
        randof <- formula(paste0(rando, '+ped(PEDIGREE_NAME)'))
      }else{
        randof <- formula('~ped(PEDIGREE_NAME)')
      }
      print(randof)
      fix <- formula('TRAIT_VALUE~1')
      
      # pd <- unique(data.frame(bdat$PEDIGREE_NAME, bdat$P1, bdat$P2))
      pd <- bdat %>% 
        select(PEDIGREE_NAME, P1, P2) %>% 
        unique() %>% 
        as.data.frame()
      pd <- pd[!duplicated(pd[,1]),]
      names(pd) <- c('Pedigree', 'P1', 'P2')
      pd$Selfing <- rep(0, rep = nrow(pd))
      
      outbred <- c('Carrot', 'Onion')
      if (crop %in% outbred){
        pd.aniv=asreml.Ainverse(pd)
      } else{
        pd.aniv=asreml.Ainverse(pd,fgen = c('Selfing',5))
      }
      
      # pd.aniv <- asreml.Ainverse(pd, fgen = c('Selfing',5))
      # pd.aniv <- asreml.Ainverse(pd)
      
      gca_mdl <- do.call(asreml, args = list(fixed = fix, random = randof, data = bdat, 
                                             na.method.X = 'include', trace = FALSE, 
                                             ginverse = list(PEDIGREE_NAME = pd.aniv$ginv),
                                             workspace = 320e+6, pworkspace = 320e+6))
      if (gca_mdl$last.message == 'LogLikelihood not converged'){
        gca_mdl <- update.asreml(gca_mdl)
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(PEDIGREE_NAME)', max.iter = 1)$predictions$pvals)
      }else{
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(PEDIGREE_NAME)', max.iter = 1)$predictions$pvals)
      }
      
      vce <- summary(gca_mdl)$varcomp
      print(vce)
      
      pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(PEDIGREE_NAME)', max.iter = 1)$predictions$pvals)
      ped_pos <- which(rownames(vce) == "ped(PEDIGREE_NAME)!ped")
      if (median(pred_df$predicted.value)==0 | (sum(c('Boundary','?','Singular') %in% vce[-ped_pos,5]) >= 1 & !vce[ped_pos,5] %in% c('Boundary','?','Singular'))) {
        gca_mdl <- rmboundary.asrtests(asrtests(gca_mdl))[[1]]
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'ped(PEDIGREE_NAME)', max.iter = 1)$predictions$pvals)
        vce <- summary(gca_mdl)$varcomp
      }
      
      h2 <- round(vce['ped(PEDIGREE_NAME)!ped',2]/sum(vce[,2]),4)
      ll <- summary(gca_mdl)$loglik
      print(c("h2: ", h2, "ll: ", ll))
      
      int_val <- as.integer(bdat$TRAIT_VALUE)
      cat_check <- sum(bdat$TRAIT_VALUE - int_val,na.rm=T)
      if (max(bdat$TRAIT_VALUE, na.rm = T) > 9){
        cat_check <- 1
      }
      
      if ((sd(pred_df$predicted.value, na.rm = T)==0)|(cat_check == 0 & (min(pred_df$predicted.value, na.rm = T)<1 | (max(bdat$TRAIT_VALUE, na.rm = T) - max(pred_df$predicted.value, na.rm = T)) > 3.5))){
        print('Poisson Distribution')
        ####CHANGING TO POISSON
        print(c('running poisson, cat_check:',cat_check,'min_pred',min(pred_df$predicted.value, na.rm = T),'max_pred:',max(pred_df$predicted.value,na.rm = T),'sd:',sd(pred_df$predicted.value, na.rm = T)))  
        gca_mdl_pois <- do.call(asreml, args = list(fixed = fix, random = randof, data = bdat, 
                                                    na.method.X = 'include', 
                                                    trace = FALSE,
                                                    family = asreml.poisson(),
                                                    ginverse = list(PEDIGREE_NAME = pd.aniv$ginv),
                                                    workspace = 320e+6, pworkspace = 320e+6))
        print(summary(gca_mdl)$varcomp)
        
        pred_df_pois <- predict(gca_mdl,classify = 'ped(PEDIGREE_NAME)',max.iter=1)$predictions$pvals
        
        if ( (sd(pred_df_pois$predicted.value, na.rm = T)>0.001) & min(pred_df_pois$predicted.value, na.rm = T) >= 1 ) {
          pred_df[,1:3] <- pred_df_pois[,1:3]
          gca_mdl <- gca_mdl_pois
        } ## IF POISS
      }
      
      rnd_ebv <- data.frame(coef(gca_mdl)$random, gca_mdl$vcoeff$random)
      names(rnd_ebv) <- c('BLUP', 'BLUP.se')
      rnd_ebv <- rnd_ebv[grep('ped', rownames(rnd_ebv)),]
      rnd_ebv$pedigree <- gsub('ped(PEDIGREE_NAME)_', '', rownames(rnd_ebv), fixed = TRUE)
      rnd_ebv$BLUP <- round(rnd_ebv$BLUP,4)
      PEV <- rnd_ebv[,2]*gca_mdl$sigma2
      # reliability <- round(1-PEV/vce['ped(PEDIGREE_NAME)!ped',2],4)
      reliability=round(1-PEV/((1 + pd.aniv$inbreeding)*vce["ped(PEDIGREE_NAME)!ped",2]),4)
      accuracy <- sqrt(reliability)
      rnd_ebv$reliability <- reliability
      
      pred_ebv <- merge(pred_df, rnd_ebv, by.x = 'PEDIGREE_NAME', by.y = 'pedigree', all.x = TRUE, all.y = TRUE)
      pred_ebv$trait <- trait
      pred_ebv$n <- NA
      pred_ebv$N <- nrow(bdat)
      pred_ebv$h2 <- h2
      pred_ebv$likelihood <- ll
      pred_ebv$mdl <- 'ped'
      n_obs <- summaryBy(TRAIT_VALUE ~ PEDIGREE_NAME, data = bdat, FUN = function(x) {length(x)})[,2]
      pred_ebv[rownames(subset(pred_ebv, as.character(PEDIGREE_NAME) %in% unique(bdat$PEDIGREE_NAME))),]$n <- n_obs
      
      # pred_ebv <- pred_ebv %>% 
      #   mutate(mdl = if_else(as.character(PEDIGREE_NAME) %in% unique(bdat$PEDIGREE_NAME[bdat$LTYPE == 'Hybrid']), 'P_SCA', 'P_GCA'))
      # 

      pred_ebv[rownames(subset(pred_ebv, as.character(PEDIGREE_NAME) %in% unique(bdat$PEDIGREE_NAME)[bdat$LTYPE == 'Hybrid'])),]$mdl <- 'P_SCA'
      pred_ebv[rownames(subset(pred_ebv, !as.character(PEDIGREE_NAME) %in% unique(bdat$PEDIGREE_NAME)[bdat$LTYPE == 'Hybrid'])),]$mdl <- 'P_GCA'
      
      freq_parent <- table(c(bdat$P1, bdat$P2))
      unique_p <- unique(c(bdat$P1, bdat$P2))
      unique_p <- unique_p[!is.na(unique_p)]
      parent_hybrid <- data.frame('Parent' = unique_p, 'Number_Crosses' = 0)
      for (u_p in 1:length(unique_p)){
        hybrid_temp <- bdat$PEDIGREE_NAME[which(bdat$P1 == unique_p[u_p] | bdat$P2 == unique_p[u_p])]
        parent_hybrid[u_p, 'Number_Crosses'] <- length(unique(hybrid_temp))
      }
      rownames(parent_hybrid) <- parent_hybrid$Parent
      pred_ebv[pred_ebv$mdl == 'P_GCA', ]$n <- parent_hybrid[as.character(pred_ebv[pred_ebv$mdl == 'P_GCA',]$PEDIGREE_NAME), 'Number_Crosses']
      if (i == 1){
        outDF <- pred_ebv
      }else{
        outDF <- rbind.data.frame(outDF, pred_ebv)
      } 
    }
  }
  GCA_result <- subset(outDF, mdl == 'P_GCA')
  SCA_result <- subset(outDF, mdl == 'P_SCA')
  blup_results <- list(GCA_result, SCA_result)
}


iblup_fun <- function(dataset, crop){
  outDF <- data.frame()
  trtList <- unique(dataset$OBSRVTN_REF_CD)
  dataset$TREP <- paste(dataset$REP_NUMBER,dataset$TRACK_ID,sep='_')
  for (i in 1:length(trtList)){
    # Some traits are not in the dataset
    # Check the existence of trait
    trait <- trtList[i]
    print(trait)
    bdat <- dataset %>% 
      filter(OBSRVTN_REF_CD == trait) %>% 
      mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
      mutate(GROWSEASON = as.character(GROWSEASON)) %>% 
      mutate(PEDIGREE_NAME = as.character(PEDIGREE_NAME)) %>% 
      mutate(BR_FIELD_ID = as.character(BR_FIELD_ID)) %>% 
      mutate(ORIGIN = as.character(ORIGIN))
    
    if (nrow(bdat) > 100){
      cols_chk <- c('BR_FIELD_ID', 'GROWSEASON', 'REPETITION', 'TREP')
      bdat[cols_chk] <- lapply(bdat[cols_chk], factor)
      levs <- sapply(bdat[cols_chk],function(x)length(unique(x)))
      print(levs)
      rando <- ""
      if (levs['BR_FIELD_ID'] > 1){
        rp <- sapply(unlist(lapply(split(bdat,bdat$BR_FIELD_ID,drop=T),function(x)length(levels(x$TREP)))),function(x) sum(x> 1))
        if(sum(rp, na.rm = T) >= 1){
          rando <- '~BR_FIELD_ID:TREP'
        }else{
          rando <- '~BR_FIELD_ID'
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
      # 
      # if (levs['TREP'] > 1){
      #   if (rando != ''){
      #     rando <- paste0(rando, '+ TREP')
      #   }else{
      #     rando <- '~ TREP'
      #   }
      # }
      if (rando != ""){
        randof <- formula(paste0(rando, '+ PEDIGREE_NAME'))
      }else{
        randof <- formula('~PEDIGREE_NAME')
      }
      print(randof)
      fix <- formula('TRAIT_VALUE~1')
      

      
      gca_mdl <- do.call(asreml, args = list(fixed = fix, random = randof, data = bdat, 
                                             na.method.X = 'include', trace = FALSE, 
                                             workspace = 320e+6, pworkspace = 320e+6))
      if (gca_mdl$last.message == 'LogLikelihood not converged'){
        gca_mdl <- update.asreml(gca_mdl)
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'PEDIGREE_NAME', max.iter = 1)$predictions$pvals)
      }else{
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'PEDIGREE_NAME', max.iter = 1)$predictions$pvals)
      }
      
      vce <- summary(gca_mdl)$varcomp
      print(vce)
      
      pred_df <- as.data.frame(predict(gca_mdl, classify = 'PEDIGREE_NAME', max.iter = 1)$predictions$pvals)
      ped_pos <- which(rownames(vce) == "PEDIGREE_NAME!ped")
      # if (median(pred_df$predicted.value)==0 | (sum(c('Boundary','?','Singular') %in% vce[-ped_pos,5]) >= 1 & !vce[ped_pos,5] %in% c('Boundary','?','Singular'))) {
        gca_mdl <- rmboundary.asrtests(asrtests(gca_mdl))[[1]]
        pred_df <- as.data.frame(predict(gca_mdl, classify = 'PEDIGREE_NAME', max.iter = 1)$predictions$pvals)
        vce <- summary(gca_mdl)$varcomp
      }
      
      int_val <- as.integer(bdat$TRAIT_VALUE)
      cat_check <- sum(bdat$TRAIT_VALUE - int_val,na.rm=T)
      if (max(bdat$TRAIT_VALUE, na.rm = T) > 9){
        cat_check <- 1
      }
      
      if ((sd(pred_df$predicted.value, na.rm = T)==0)|(cat_check == 0 & (min(pred_df$predicted.value, na.rm = T)<1 | (max(bdat$TRAIT_VALUE, na.rm = T) - max(pred_df$predicted.value, na.rm = T)) > 3.5))){
        print('Poisson Distribution')
        ####CHANGING TO POISSON
        print(c('running poisson, cat_check:',cat_check,'min_pred',min(pred_df$predicted.value, na.rm = T),'max_pred:',max(pred_df$predicted.value,na.rm = T),'sd:',sd(pred_df$predicted.value, na.rm = T)))  
        gca_mdl_pois <- do.call(asreml, args = list(fixed = fix, random = randof, data = bdat, 
                                                    na.method.X = 'include', 
                                                    trace = FALSE,
                                                    family = asreml.poisson(),
                                                    workspace = 320e+6, pworkspace = 320e+6))
        print(summary(gca_mdl)$varcomp)
        
        pred_df_pois <- predict(gca_mdl,classify = 'PEDIGREE_NAME',max.iter=1)$predictions$pvals
        
        if ( (sd(pred_df_pois$predicted.value, na.rm = T)>0.001) & min(pred_df_pois$predicted.value, na.rm = T) >= 1 ) {
          pred_df[,1:3] <- pred_df_pois[,1:3]
          gca_mdl <- gca_mdl_pois
        } ## IF POISS
      }
      
      pred_df$trait <- trait
      
      if (i == 1){
        outDF <- pred_df 
      }else{
        outDF <- rbind(outDF, pred_df)
      }
  }
  return(outDF)
  }

