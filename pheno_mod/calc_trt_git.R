calc_trt=function(pheno,Crop,prog){
  crop=Crop
  # pheno$TRAIT_VALUE=as.numeric(pheno$TRAIT_VALUE)
  if (crop == 'Onion'){
    if (prog %in% c('RV', 'SL', 'UH')){
      pheno$subplot_id <- paste(pheno$PLOT_BID, pheno$REPETITION, sep = "_")
      if (sum(c('TBLNO', 'STAND') %in% unique(pheno$OBSRVTN_REF_CD)) == 2){  
        pheno_dt=subset(pheno,OBSRVTN_REF_CD=='TBLNO')
        pheno_ct=subset(pheno,OBSRVTN_REF_CD !='TBLNO')
        pheno_ds=subset(pheno,OBSRVTN_REF_CD=='STAND')
        # pheno_ds=subset(pheno_ds,PLOT_ID %in% pheno_dt$PLOT_ID)
        colnames(pheno_ds)[which(colnames(pheno_ds) == 'TRAIT_VALUE')] <- 'STAND'
        pheno_combine <- pheno_dt %>% inner_join(pheno_ds[,c('subplot_id', 'STAND')], by = 'subplot_id')
        dim(pheno_combine)
        # dim(pheno_ds)
        if (nrow(pheno_combine) > 0){
          stnd=as.numeric(pheno_combine$STAND)
          blt=as.numeric(pheno_combine$TRAIT_VALUE)
          pheno_tc=pheno_combine[,-ncol(pheno_combine)]
          pheno_tc$OBSRVTN_REF_CD='TBLNO_C'
          pheno_tc$FTS_OBSERVATION="Bolting corrected by Stand"
          pheno_tc$TRAIT_VALUE=NA
          pheno_tc$TRAIT_VALUE=blt/stnd
          pheno_dt=rbind(pheno_dt,pheno_tc)
          pheno=rbind(pheno_ct,pheno_dt) 
        }
      }
      # if (sum(c('YIELDWT', 'STAND', 'HAARE') %in% unique(pheno$OBSRVTN_REF_CD)) == 3){
      if (sum(c('STAND', 'HAARE') %in% unique(pheno$OBSRVTN_REF_CD)) == 2 & length(grep('YIELDWT', unique(pheno$OBSRVTN_REF_CD)))>=1){
        pheno_d <- pheno[grep('YIELDWT', pheno$OBSRVTN_REF_CD),]
        pheno_c <- pheno[!grepl('YIELDWT', pheno$OBSRVTN_REF_CD),]
        pheno_ds=subset(pheno,OBSRVTN_REF_CD=='STAND')
        pheno_dh=subset(pheno,OBSRVTN_REF_CD=='HAARE')
        colnames(pheno_ds)[which(colnames(pheno_ds) == 'TRAIT_VALUE')] <- "STAND"
        colnames(pheno_dh)[which(colnames(pheno_dh) == 'TRAIT_VALUE')] <- "HAARE"
        pheno_cy <- pheno_d %>% left_join(pheno_ds[,c('subplot_id', 'STAND')], by = 'subplot_id') %>% 
          left_join(pheno_dh[, c('subplot_id', 'HAARE')], by = 'subplot_id')
        yld=as.numeric(pheno_cy$TRAIT_VALUE)
        stnd=as.numeric(pheno_cy$STAND)
        har=as.numeric(pheno_cy$HAARE)
        yld_har <- yld/har
        stnd_har <- stnd/har
        
        pheno_di=pheno_cy[,-which(colnames(pheno_cy) %in% c('STAND', 'HAARE'))]
        pheno_di$OBSRVTN_REF_CD='YIELDWT_HAARE'
        pheno_di$FTS_OBSERVATION="Yield corrected by area"
        pheno_di$TRAIT_VALUE=NA
        pheno_di$TRAIT_VALUE <- yld_har
        pheno_di <- pheno_di %>% group_by(subplot_id) %>% 
          mutate(TRAIT_VALUE = sum(TRAIT_VALUE, na.rm = TRUE))
        pheno_di <- as.data.frame(pheno_di[!duplicated(pheno_di[,c('subplot_id', 'PEDIGREE_NAME')]),])
        # mkyld=c('SEW2-3','SEW3-3.75','SEW>3.75',"NWW2.25-3","NWW3-4","NWW4-4.25","NWW>4.25")
        pheno_std <- pheno_cy[,-which(colnames(pheno_cy) %in% c('STAND', 'HAARE'))]
        pheno_std$OBSRVTN_REF_CD='STAND_HAARE'
        pheno_std$FTS_OBSERVATION="Stand corrected by area"
        pheno_std$TRAIT_VALUE=NA
        pheno_std$TRAIT_VALUE <- stnd_har
        pheno_std <- pheno_std %>% 
          group_by(subplot_id) %>% 
          mutate(TRAIT_VALUE = sum(TRAIT_VALUE, na.rm = TRUE))
        pheno_std <- as.data.frame(pheno_std[!duplicated(pheno_std[, c('subplot_id', 'PEDIGREE_NAME')]),])
        pheno <- do.call(rbind, list(pheno_c, pheno_d,pheno_di, pheno_std))
      }
    }
  }
  
  if (crop == 'Cucumber'){ 
    pheno$TRAIT_VALUE=as.numeric(pheno$TRAIT_VALUE)
    if (prog == 'Y3'){
      if ('NETWT' %in% unique(pheno$OBSRVTN_REF_CD)){
        posi=which(names(pheno)=='TRAIT_VALUE')
        #Netwt
        pheno_c_nt=subset(pheno,OBSRVTN_REF_CD=='NETWT')
        pheno_c_nt$ntrt=interaction(pheno_c_nt$PLOT_ID,pheno_c_nt$GERMPLASM_ID,pheno_c_nt$REPETITION)
        names(pheno_c_nt)[posi]='NT'
        #MKYld
        pheno_c_fk=subset(pheno,OBSRVTN_REF_CD=='FRNMK')
        pheno_c_fk$ntrt=interaction(pheno_c_fk$PLOT_ID,pheno_c_fk$GERMPLASM_ID,pheno_c_fk$REPETITION)
        names(pheno_c_fk)[posi]='FK'
        #NPLT
        pheno_c_np=subset(pheno,OBSRVTN_REF_CD=='NHPLT')##Not the same size. Take Plant number more than yield measures.
        pheno_c_np$ntrt=interaction(pheno_c_np$PLOT_ID,pheno_c_np$GERMPLASM_ID,pheno_c_np$REPETITION)
        names(pheno_c_np)[posi]='NP'
        #NPA (Density)
        pheno_c_py=subset(pheno,OBSRVTN_REF_CD=='NPA')  #USED TO BE YIELD NPA
        pheno_c_py$ntrt=interaction(pheno_c_py$PLOT_ID,pheno_c_py$GERMPLASM_ID,pheno_c_py$REPETITION)
        names(pheno_c_py)[posi]='PY'
        dim(pheno_c_py)
        #Combine data based on plot,gpid,trep + rep
        pheno_comb=Reduce(function(x, y) merge(x, y,by='ntrt'), list(pheno_c_nt,pheno_c_fk,pheno_c_np,pheno_c_py))
        nt_wt=as.numeric(pheno_comb$NT)
        ft_mk=as.numeric(pheno_comb$FK)
        n_pt=as.numeric(pheno_comb$NP)
        n_pa=as.numeric(pheno_comb$PY)
        #Calc traits
        #FW=kg/frt num
        # AFW_C=nt_wt/ft_mk*1000
        #FRNM=FT/#plt*density
        NFSA_C=ft_mk/n_pt*n_pa
        #NtWt =KG/Plt#*Density
        WFSA_C=nt_wt/n_pt*n_pa
        pheno_comb=pheno_comb[,2:(dim(pheno)[2]+1)]
        names(pheno_comb)=names(pheno)
        # pheno_c_nt=pheno_comb
        # pheno_c_nt$OBSRVTN_REF_CD='AFW_C'
        # pheno_c_nt$TRAIT_VALUE=AFW_C
        
        nt_dat <- pheno_c_nt %>% 
          dplyr::group_by(PLOT_BID) %>% 
          dplyr::summarise(sum_nt = sum(NT, rm.na = T))
        
        
        fk_dat <- pheno_c_fk %>% 
          dplyr::group_by(PLOT_BID) %>% 
          dplyr::summarise(sum_fk = sum(FK, rm.na = T))
        
        nt_fk <- nt_dat %>% 
          inner_join(fk_dat, by = 'PLOT_BID') %>% 
          mutate(AFW_C = sum_nt/sum_fk)
        
        pheno_c_nt <- pheno_comb[!duplicated(pheno_comb$PLOT_BID),]
        
        pheno_c_nt <- pheno_c_nt %>% 
          inner_join(nt_fk, by = 'PLOT_BID') %>% 
          mutate(TRAIT_VALUE = AFW_C,
                 OBSRVTN_REF_CD = 'AFW_C',
                 REPETITION = 1) %>% 
          dplyr::select(-c('sum_nt', 'sum_fk', 'AFW_C'))
        
        pheno_c_fk=pheno_comb
        pheno_c_fk$OBSRVTN_REF_CD='NFSA_C'
        pheno_c_fk$TRAIT_VALUE=NFSA_C
        pheno_c_py=pheno_comb
        pheno_c_py$OBSRVTN_REF_CD='WFSA_C'
        pheno_c_py$TRAIT_VALUE=WFSA_C
        pheno = rbind(pheno,pheno_c_nt,pheno_c_fk,pheno_c_py)
      }
    }
    
    ### Calculated trait Nischit program
    if (prog %in% c('M3','EF')){
      pheno$subplot_id <- paste(pheno$PLOT_BID, pheno$REPETITION, sep = "_")
      pheno$TRAIT_VALUE<-as.numeric(as.character(pheno$TRAIT_VALUE))
      
      #if (sum(c('TBLNO', 'STAND') %in% unique(pheno$OBSRVTN_REF_CD)) == 2){  check col names
      if (sum(c("YLD1A" , "YLD1B", "YLD2A","YLD2B" ,"YLD3A" , "YLD3B", "YLD4" , "YLDC","FDAYHV","STAND") %in% unique(pheno$OBSRVTN_REF_CD)) == 10){  #check col names
        
        pheno_1a=subset(pheno,OBSRVTN_REF_CD=='YLD1A')
        pheno_1b=subset(pheno,OBSRVTN_REF_CD=='YLD1B')
        pheno_2a=subset(pheno,OBSRVTN_REF_CD=='YLD2A')
        pheno_2b=subset(pheno,OBSRVTN_REF_CD=='YLD2B')
        pheno_3a=subset(pheno,OBSRVTN_REF_CD=='YLD3A')
        pheno_3b=subset(pheno,OBSRVTN_REF_CD=='YLD3B')
        pheno_4y=subset(pheno,OBSRVTN_REF_CD=='YLD4')
        pheno_c=subset(pheno,OBSRVTN_REF_CD=='YLDC')
        pheno_fday=subset(pheno,OBSRVTN_REF_CD=='FDAYHV')
        pheno_std=subset(pheno,OBSRVTN_REF_CD=='STAND')
        
        colnames(pheno_1a)[which(colnames(pheno_1a) == 'TRAIT_VALUE')] <- 'YLD1A'
        colnames(pheno_1b)[which(colnames(pheno_1b) == 'TRAIT_VALUE')] <- 'YLD1B'
        colnames(pheno_2a)[which(colnames(pheno_2a) == 'TRAIT_VALUE')] <- 'YLD2A'
        colnames(pheno_2b)[which(colnames(pheno_2b) == 'TRAIT_VALUE')] <- 'YLD2B'
        colnames(pheno_3a)[which(colnames(pheno_3a) == 'TRAIT_VALUE')] <- 'YLD3A'
        colnames(pheno_3b)[which(colnames(pheno_3b) == 'TRAIT_VALUE')] <- 'YLD3B'
        colnames(pheno_4y)[which(colnames(pheno_4y) == 'TRAIT_VALUE')] <- 'YLD4'
        colnames(pheno_c)[which(colnames(pheno_c) == 'TRAIT_VALUE')] <- 'YLDC'
        colnames(pheno_fday)[which(colnames(pheno_fday) == 'TRAIT_VALUE')] <- 'FDAYHV'
        colnames(pheno_std)[which(colnames(pheno_std) == 'TRAIT_VALUE')] <- 'STAND'
        
        pheno_combine <- pheno_1a %>% inner_join(pheno_1b[,c('subplot_id', 'YLD1B')], by = 'subplot_id')
        pheno_combine <- pheno_combine %>% inner_join(pheno_2a[,c('subplot_id', 'YLD2A')], by = 'subplot_id')
        pheno_combine <- pheno_combine %>% inner_join(pheno_2b[,c('subplot_id', 'YLD2B')], by = 'subplot_id')
        pheno_combine <- pheno_combine %>% inner_join(pheno_3a[,c('subplot_id', 'YLD3A')], by = 'subplot_id')
        pheno_combine <- pheno_combine %>% inner_join(pheno_3b[,c('subplot_id', 'YLD3B')], by = 'subplot_id')
        pheno_combine <- pheno_combine %>% inner_join(pheno_4y[,c('subplot_id', 'YLD4')], by = 'subplot_id')
        pheno_combine <- pheno_combine %>% inner_join(pheno_c[,c('subplot_id', 'YLDC')], by = 'subplot_id')
        pheno_combine <- pheno_combine %>% inner_join(pheno_fday[,c('subplot_id', 'FDAYHV')], by = 'subplot_id')
        pheno_combine <- pheno_combine %>% inner_join(pheno_std[,c('subplot_id', 'STAND')], by = 'subplot_id')
        
        # dim(pheno_ds)
        if (nrow( pheno_combine) > 0){
          #tryCatch({ 
          pheno_combine <- pheno_combine%>%rowwise()%>%mutate(total=sum(YLD1A,YLD1B,YLD2A,YLD2B,YLD3A,YLD3B,YLD4,YLDC))
          pheno_combine<-pheno_combine%>%mutate(fiveper_C=YLD4-total*5/100) 
          #pheno_dt = pheno[!pheno$OBSRVTN_REF_CD%in%c("YLD1A","YLD1B","YLD2A","YLD2B","YLD3A","YLD3B","YLD4","YLDC","YLD2"),]
          pheno_combine = data.frame(pheno_combine)
          pheno_combine=  droplevels(pheno_combine,pheno_combine$SUB_COUNTRY)
          
          location = unique(pheno_combine$SUB_COUNTRY)
          weight_final = data.frame()
          
          for (i in 1:length(location)){
            #i=1
            pheno_combine_1<-pheno_combine[pheno_combine$SUB_COUNTRY==location[[i]],]
            pheno_combine_1 = pheno_combine_1%>%droplevels("SUB_COUNTRY")
            pheno_combine_2= pheno_combine%>%mutate(CROSS_REP=paste0(CROSS_NAME,TREP))
            
            mt_fit <- suppressWarnings(asreml(cbind(YLD1A,YLD1B,YLD2A,YLD2B,YLD3A,YLD3B)~trait,
                                              random =~us(trait):CROSS_REP+fiveper_C + YLDC  + STAND + FDAYHV,
                                              na.method.X='include',trace=FALSE,family=asreml.gaussian(),rcov=~units:us(trait),
                                              data = pheno_combine_2,
                                              workspace = 320e+6, pworkspace = 320e+6))
            
            results<-data.frame(summary(mt_fit,all=TRUE)$coef.random)
            res<-data.frame("ID"=rownames(results),results,row.names = NULL)
            res$ID<-gsub("trait_","",res$ID)
            #res$ID<-gsub("CROSS_REP_","",res$ID)
            res<-res[grep("YLD", res$ID),]
            res<-res[!res$ID%in%c("YLDC"),]
            res<-res%>%separate(ID,c('TRAIT','CROSS_NAME'),sep=':') 
            res<-res[!colnames(res)%in%c("std.error","z.ratio")]
            res1<-res%>%group_by(CROSS_NAME)%>%spread(TRAIT,solution)
            res2<-data.frame("CROSS_REP"=res1$CROSS_NAME,"YLD1A_C"=res1$YLD1A+mean(pheno_combine_2$YLD1A),"YLD1B_C"=res1$YLD1B+mean(pheno_combine_2$YLD1B),
                             "YLD2A_C"=res1$YLD2A+mean(pheno_combine_2$YLD2A),"YLD2B_C"=res1$YLD2B+mean(pheno_combine_2$YLD2B),
                             "YLD3A_C"=res1$YLD3A+mean(pheno_combine_2$YLD3A),"YLD3B_C"=res1$YLD3B+mean(pheno_combine_2$YLD3B)) # prediction in weight scale 
            res2$YLD1_3_C = res2$YLD1A_C+res2$YLD1B_C+res2$YLD2A_C+res2$YLD2B_C+res2$YLD3A_C+res2$YLD3B_C
            res2$CROSS_REP<-gsub("CROSS_REP_","",res2$CROSS_REP)
            res_weight_loc<-suppressMessages(left_join(pheno_combine_2,res2))
            weight_final<-rbind(weight_final,res_weight_loc)
          }
          
          pheno_tc= weight_final[,!colnames(weight_final)%in%c("YLD1A","YLD1B" ,"YLD2A" ,"YLD2B", "YLD3A" , "YLD3B" ,"YLD4","YLDC", "FDAYHV" , "STAND","total" , "fiveper_C" , "paste0(CROSS_NAME, REPETITION)", "CROSS_REP","YLD1A_C",
                                                               "YLD1B_C", "YLD2A_C", "YLD2B_C" , "YLD3A_C" ,"YLD3B_C")]
          pheno_tc$OBSRVTN_REF_CD='YLD1_3_C'
          pheno_tc$FTS_OBSERVATION="total weight sum grades 1- 3"
          pheno_tc$TRAIT_VALUE=weight_final$YLD1_3_C
          pheno_tc = pheno_tc[,!colnames(pheno_tc)%in%c("YLD1_3_C")]
          pheno=rbind(pheno,pheno_tc) 
          pheno = pheno [,-ncol(pheno)]
        }   
      }
    }
  } 
  
  # Calcualte trait for squash. Program is R5
  if (crop == "Squash"){
    if (prog == "R5"){
      pheno$TRAIT_VALUE=as.numeric(pheno$TRAIT_VALUE)
      unique_Ob <- unique(pheno$OBSRVTN_REF_CD)
      posi <- which(names(pheno)=='TRAIT_VALUE')
      pheno$OBS_DATE_RECORDED <- as.Date(as.character(pheno$OBS_DATE_RECORDED), '%Y-%m-%d')
      # pheno$OBS_DATE_RECORDED <- as.Date(as.character(pheno$OBS_DATE_RECORDED))
      
      # NETWT: weight of marketable fresh fruit
      # FRNMK: number of marketable fruit harvested
      if (length(grep('NETWT', unique_Ob))>=1){
        if (length(grep('FRNMK', unique_Ob))>=1){
          pheno_c_netwt <- pheno[grep('NETWT', pheno$OBSRVTN_REF_CD), ]
          pheno_c_netwt$ntrt <- interaction(pheno_c_netwt$PLOT_ID,pheno_c_netwt$GERMPLASM_ID,pheno_c_netwt$REPETITION)
          names(pheno_c_netwt)[posi] <- 'NT'
          
          pheno_c_frnmk <- pheno[grep('FRNMK', pheno$OBSRVTN_REF_CD), ]
          pheno_c_frnmk$ntrt <- interaction(pheno_c_frnmk$PLOT_ID,pheno_c_frnmk$GERMPLASM_ID,pheno_c_frnmk$REPETITION)
          names(pheno_c_frnmk)[posi] <- 'FR'
          
          # Calculate Trait
          if (length(intersect(pheno_c_netwt$ntrt, pheno_c_frnmk$ntrt))>=1){
            pheno_nt_fr <- merge(pheno_c_netwt, pheno_c_frnmk, by = 'ntrt')
            n_netwt <- as.numeric(pheno_nt_fr$NT)
            n_frnmk <- as.numeric(pheno_nt_fr$FR)
            
            # Average weight per marketable fresh fruit harvested
            AFW_H <- n_netwt/n_frnmk
            
            pheno_c_netwt <- pheno_nt_fr[,2:(dim(pheno)[2]+1)]
            names(pheno_c_netwt) <- names(pheno)
            pheno_c_netwt$OBSRVTN_REF_CD <- 'AFW_H'
            pheno_c_netwt$TRAIT_VALUE <- AFW_H
            
            pheno <- rbind(pheno, pheno_c_netwt)
          }
        }
      }
      
      
      
      
      
      # YNMPL: number of marketable fruit per plant
      # WMFPL: weight of marketable fruit per plant
      if (length(grep('YNMPL', unique_Ob))>=1){
        if (length(grep('WMFPL', unique_Ob))>=1){
          pheno_c_ynmpl <- pheno[grep('YNMPL', pheno$OBSRVTN_REF_CD), ]
          pheno_c_ynmpl$ntrt <- interaction(pheno_c_ynmpl$PLOT_ID,pheno_c_ynmpl$GERMPLASM_ID,pheno_c_ynmpl$REPETITION)
          names(pheno_c_ynmpl)[posi] <- 'YL'
          
          pheno_c_wmfpl <- pheno[grep('WMFPL', pheno$OBSRVTN_REF_CD), ]
          pheno_c_wmfpl$ntrt <- interaction(pheno_c_wmfpl$PLOT_ID,pheno_c_wmfpl$GERMPLASM_ID,pheno_c_wmfpl$REPETITION)
          names(pheno_c_wmfpl)[posi] <- 'WL'
          
          # Calculate Trait
          if (length(intersect(pheno_c_wmfpl$ntrt, pheno_c_ynmpl$ntrt))>=1){
            pheno_wl_yl <- merge(pheno_c_wmfpl, pheno_c_ynmpl, by = "ntrt")
            n_ynmpl <- as.numeric(pheno_wl_yl$YL)
            n_wmfpl <- as.numeric(pheno_wl_yl$WL)
            # Average weight per marketable fruit per plant
            AFW_PL <- n_wmfpl/n_ynmpl
            
            pheno_c_wmfpl <- pheno_wl_yl[,2:(dim(pheno)[2]+1)]
            names(pheno_c_wmfpl) <- names(pheno)
            pheno_c_wmfpl$OBSRVTN_REF_CD <- 'AFW_PL'
            pheno_c_wmfpl$TRAIT_VALUE <- AFW_PL
            
            pheno <- rbind(pheno, pheno_c_wmfpl)
          }
        }
      }
      
      # WMFSA: weight of marketable fresh furit per square meter/foot
      # NPA: number of plants per square meter/foot
      # NHPLT: number of harvested plants per plot
      if (length(grep('WMFSA', unique_Ob))>=1){
        if (length(grep('NPA', unique_Ob))>=1){
          if (length(grep('NHPLT', unique_Ob))>=1){
            pheno_c_wmfsa <- pheno[grep('WMFSA', pheno$OBSRVTN_REF_CD), ]
            pheno_c_wmfsa$ntrt <- interaction(pheno_c_wmfsa$PLOT_ID,pheno_c_wmfsa$GERMPLASM_ID,pheno_c_wmfsa$REPETITION)
            names(pheno_c_wmfsa)[posi] <- 'WA'
            
            pheno_c_npa <- pheno[grep('NPA', pheno$OBSRVTN_REF_CD), ]
            pheno_c_npa$ntrt <- interaction(pheno_c_npa$PLOT_ID,pheno_c_npa$GERMPLASM_ID,pheno_c_npa$REPETITION)
            names(pheno_c_npa)[posi] <- 'NP'
            
            pheno_c_nhplt <- pheno[grep('NHPLT', pheno$OBSRVTN_REF_CD), ]
            pheno_c_nhplt$ntrt <- interaction(pheno_c_nhplt$PLOT_ID, pheno_c_nhplt$GERMPLASM_ID, pheno_c_nhplt$REPETITION)
            names(pheno_c_nhplt)[posi] <- 'NH' 
            
            # Calculate trait
            if (length(Reduce(f = function(x,y) intersect(x,y), 
                              list(pheno_c_wmfsa, pheno_c_npa, pheno_c_nhplt)))>=1){
              pheno_wa_np_nh <- Reduce(f = function(x, y) merge(x, y, by  = 'ntrt'), 
                                       list(pheno_c_wmfsa, pheno_c_npa, pheno_c_nhplt))
              n_wmfsa <- as.numeric(pheno_wa_np_nh$WA)
              n_npa <- as.numeric(pheno_wa_np_nh$NP)
              n_nhplt <- as.numeric(pheno_wa_np_nh$NH)
              
              # Average weight per fruit per square meter/foot
              AFW_NPA <- n_wmfsa/((n_wmfsa/n_nhplt)*n_npa)
              
              pheno_c_wmfsa <- pheno_wa_np_nh[,2:(dim(pheno)[2]+1)]
              names(pheno_c_wmfsa) <- names(pheno)
              pheno_c_wmfsa$OBSRVTN_REF_CD <- 'AFW_NPA'
              pheno_c_wmfsa$TRAIT_VALUE <- AFW_NPA
              
              pheno <- rbind(pheno, pheno_c_wmfsa)
            }
          }
        }
      }
      
      # WMFFP: weight of marketable fresh fruit per plot
      # YNMPP: number of marketable fruits per plot
      # NHPLT: number of haversted plants per plot
      if (length(grep('NHPLT', unique_Ob))>=1){
        pheno_c_nhplt <- pheno[grep('NHPLT', pheno$OBSRVTN_REF_CD), ]
        # Check if NHPLT changes over the trial
        chk_nhplt <- aggregate(OBSRVTN_REF_CD~PLOT_ID + REPETITION + GERMPLASM_ID, data = pheno_c_nhplt, FUN = function(x){length(unique(x))})
        if (sum(chk_nhplt[2] == 1) == dim(chk_nhplt)[1]){
          # Case 1
          if (length(grep('YNMPP', unique_Ob))>=1){
            if (length(grep('WMFFP', unique_Ob))>=1){
              pheno_c_nhplt$ntrt <- interaction(pheno_c_nhplt$PLOT_ID, pheno_c_nhplt$GERMPLASM_ID)
              names(pheno_c_nhplt)[posi] <- 'NH'
              
              pheno_c_ynmpp <- pheno[grep('YNMPP', pheno$OBSRVTN_REF_CD), ]
              pheno_c_ynmpp$ntrt <- interaction(pheno_c_ynmpp$PLOT_ID, pheno_c_ynmpp$GERMPLASM_ID)
              names(pheno_c_ynmpp)[posi] <- 'YP'
              
              pheno_c_wmffp <- pheno[grep('WMFFP', pheno$OBSRVTN_REF_CD), ]
              pheno_c_wmffp$ntrt <- interaction(pheno_c_wmffp$PLOT_ID, pheno_c_wmffp$GERMPLASM_ID)
              names(pheno_c_wmffp)[posi] <- 'WP'
              pheno_c_ynmpp$YP <- as.numeric(pheno_c_ynmpp$YP)
              # pheno_afw_pp$WP <- as.numeric(pheno_afw_pp$WP)
              # YNMPL
              pheno_yp_nl <- merge(pheno_c_ynmpp,pheno_c_nhplt, by = "ntrt")
              n_ynmpp <- as.numeric(pheno_yp_nl$YP)
              n_nhplt <- as.numeric(pheno_yp_nl$NH)
              
              # Number of fruit per plant
              YNMPL <- n_ynmpp/n_nhplt
              
              pheno_yp_nl <- pheno_yp_nl[,1:dim(pheno)[2]+1]
              names(pheno_yp_nl) <- names(pheno)
              pheno_yp_nl$OBSRVTN_REF_CD <- 'YNMPL'
              pheno_yp_nl$TRAIT_VALUE <- YNMPL
              
              pheno <- rbind(pheno, pheno_yp_nl)
              # WMFPL
              pheno_c_ynmpp$ntrt_1 <- interaction(pheno_c_ynmpp$PLOT_ID, pheno_c_ynmpp$GERMPLASM_ID, pheno_c_ynmpp$OBS_DATE_RECORDED)
              pheno_c_wmffp$ntrt_1 <- interaction(pheno_c_wmffp$PLOT_ID, pheno_c_wmffp$GERMPLASM_ID, pheno_c_wmffp$OBS_DATE_RECORDED)
              pheno_yp_wp <- merge(pheno_c_ynmpp, pheno_c_wmffp, by = 'ntrt_1')
              n_ynmpp <- aggregate(YP~ntrt.x, data = pheno_yp_wp, sum)
              n_wmffp <- as.numeric(pheno_yp_wp$WP)
              
              pheno_yp_wp_1 <- merge(pheno_yp_wp, n_ynmpp, by= 'ntrt.x', all.x = T)
              
              fruit_weight <- pheno_yp_wp_1$WP/pheno_yp_wp_1$YP.y # AFW_PL:average fruit per fruit per plant
              
              pheno_yp_wp_1 <- pheno_yp_wp_1[, 1:dim(pheno)[2]+2]
              names(pheno_yp_wp_1) <- names(pheno)
              pheno_yp_wp_1$OBSRVTN_REF_CD <- 'AFW_PL'
              pheno_yp_wp_1$TRAIT_VALUE <- fruit_weight
              
              pheno <- rbind(pheno, pheno_yp_wp_1)
              
              pheno_yp_nl$ntrt <- interaction(pheno_yp_nl$PLOT_ID, pheno_yp_nl$GERMPLASM_ID)
              ynmpl_sum <- aggregate(TRAIT_VALUE~ntrt, data = pheno_yp_nl, sum)
              pheno_yp_wp_1$ntrt <- interaction(pheno_yp_wp_1$PLOT_ID, pheno_yp_wp_1$GERMPLASM_ID)
              pheno_wp_nl <- merge(pheno_yp_wp_1, ynmpl_sum, by = 'ntrt', all.x = T)
              
              WMFPL <- pheno_wp_nl$TRAIT_VALUE.x*pheno_wp_nl$TRAIT_VALUE.y
              
              pheno_wp_nl <- pheno_wp_nl[,2:(dim(pheno_wp_nl)[2]-1)]
              names(pheno_wp_nl) <- names(pheno)
              pheno_wp_nl$OBSRVTN_REF_CD <- 'WMFPL'
              pheno_wp_nl$TRAIT_VALUE <- WMFPL
              pheno <- rbind(pheno, pheno_wp_nl)
              
              # AFW_PP: average fruit weight per fruit per plot
              ynmpp_sum <- aggregate(YP~ntrt, data = pheno_c_ynmpp, sum)
              pheno_afw_pp <- merge(pheno_c_wmffp, ynmpp_sum, by = 'ntrt', all.x = T)
              AFW_PP <- pheno_afw_pp$WP/pheno_afw_pp$YP
              pheno_afw_pp <- pheno_afw_pp[, 2:(dim(pheno_afw_pp)[2]-2)]
              names(pheno_afw_pp) <- names(pheno)
              pheno_afw_pp$OBSRVTN_REF_CD <- 'AFW_PP'
              pheno_afw_pp$TRAIT_VALUE <- AFW_PP
              pheno <- rbind(pheno, pheno_afw_pp)
              
            }
          }
        }else{
          # Case 2
          # pheno_c_nhplt$ntrt <- interaction(pheno_c_nhplt$PLOT_ID, pheno_c_nhplt$GERMPLASM_ID, 
          #                                   pheno_c_nhplt$REPETITION)
          names(pheno_c_nhplt)[posi] <- 'NH'
          dates <- sort(unique(pheno_c_nhplt$OBS_DATE_RECORDED))
          dates_interval <- list()
          for (i in 1:(length(dates))){
            if (i != length(dates)){
              seq_dates <- seq.Date(dates[i], dates[i+1]-1, by = 1)
              dates_interval[length(dates_interval)+1] <- list(seq_dates)
            }else{
              seq_dates <- seq.Date(dates[i], dates[i]+30, by = 1)
              dates_interval[length(dates_interval)+1] <- list(seq_dates)
            }
          }
          
          if (length(grep('WMFFP', unique_Ob))>=1){
            if (length(grep('YNMPP', unique_Ob))>=1){
              
              pheno_c_wmffp <- pheno[grep('WMFFP', pheno$OBSRVTN_REF_CD), ]
              pheno_c_wmffp$coded_date <- NA
              # Recode Observation Date
              for (i in 1:length(dates_interval)){
                pheno_c_wmffp$coded_date[pheno_c_wmffp$OBS_DATE_RECORDED %in% dates_interval[[i]]] <- i
                pheno_c_nhplt$coded_date[pheno_c_nhplt$OBS_DATE_RECORDED %in% dates_interval[[i]]] <- i
              }
              
              pheno_c_wmffp$ntrt <- interaction(pheno_c_wmffp$PLOT_ID, pheno_c_wmffp$GERMPLASM_ID, 
                                                pheno_c_wmffp$coded_date)
              pheno_c_nhplt$ntrt <- interaction(pheno_c_nhplt$PLOT_ID, pheno_c_nhplt$GERMPLASM_ID,
                                                pheno_c_nhplt$coded_date)
              # pheno_c_wmffp$ntrt <- interaction(pheno_c_wmffp$PLOT_ID, pheno_c_wmffp$GERMPLASM_ID, pheno_c_wmffp$REPETITION)
              names(pheno_c_wmffp)[posi] <- "WP"
              
              pheno_c_ynmpp <- pheno[grep('YNMPP', pheno$OBSRVTN_REF_CD), ]
              pheno_c_ynmpp$coded_date <- NA
              # Recode Observation Date
              for (i in 1:length(dates_interval)){
                pheno_c_ynmpp$coded_date[pheno_c_ynmpp$OBS_DATE_RECORDED %in% dates_interval[[i]]] <- i
              }
              pheno_c_ynmpp$ntrt <- interaction(pheno_c_ynmpp$PLOT_ID, pheno_c_ynmpp$GERMPLASM_ID, pheno_c_ynmpp$coded_date)
              # pheno_c_ynmpp$ntrt <- interaction(pheno_c_ynmpp$PLOT_ID, pheno_c_ynmpp$GERMPLASM_ID, pheno_c_ynmpp$REPETITION)
              names(pheno_c_ynmpp)[posi] <- 'YP'
              
              # Calcuated Trait: AFW_PP = WMFFP/YNMPP
              if (length(intersect(pheno_c_ynmpp$ntrt, pheno_c_wmffp$ntrt))>=1){
                pheno_c_wmffp$ntrt_1 <- interaction(pheno_c_wmffp$PLOT_ID, pheno_c_wmffp$GERMPLASM_ID)
                pheno_c_ynmpp$ntrt_1 <- interaction(pheno_c_ynmpp$PLOT_ID, pheno_c_ynmpp$GERMPLASM_ID)
                ynmpp_sum <- aggregate(YP~ntrt_1, data = pheno_c_ynmpp, sum)
                pheno_afw_pp <- merge(pheno_c_wmffp, ynmpp_sum, by = 'ntrt_1', all.x = T)
                AFW_PP <- pheno_afw_pp$WP/pheno_afw_pp$YP
                pheno_afw_pp <- pheno_afw_pp[, 2:(dim(pheno_afw_pp)[2]-3)]
                names(pheno_afw_pp) <- names(pheno)
                pheno_afw_pp$OBSRVTN_REF_CD <- 'AFW_PP'
                pheno_afw_pp$TRAIT_VALUE <- AFW_PP
                pheno <- rbind(pheno, pheno_afw_pp)
              }
              # Calculated Trait: YNMPL
              
              pheno_yp_nl <- pheno_c_ynmpp %>% left_join(pheno_c_nhplt, by = "ntrt")
              n_ynmpp <- as.numeric(pheno_yp_nl$YP)
              n_nhplt <- as.numeric(pheno_yp_nl$NH)
              
              # Number of fruit per plant
              YNMPL <- n_ynmpp/n_nhplt
              
              pheno_yp_nl <- pheno_yp_nl[,1:dim(pheno)[2]]
              names(pheno_yp_nl) <- names(pheno)
              pheno_yp_nl$OBSRVTN_REF_CD <- 'YNMPL'
              pheno_yp_nl$TRAIT_VALUE <- YNMPL
              
              pheno <- rbind(pheno, pheno_yp_nl)
              
              # Calculated Trait: WMFPL
              # pheno_wp_nl <- pheno_c_wmffp %>% left_join(pheno_c_nhplt, by = "ntrt")
              # n_wmffp <- as.numeric(pheno_wp_nl$WP)
              # n_nhplt <- as.numeric(pheno_wp_nl$NH)
              for (i in 1:length(dates_interval)){
                pheno_yp_nl$coded_date[pheno_yp_nl$OBS_DATE_RECORDED %in% dates_interval[[i]]] <- i
              }
              pheno_yp_nl$ntrt <- interaction(pheno_yp_nl$PLOT_ID, 
                                              pheno_yp_nl$GERMPLASM_ID, 
                                              pheno_yp_nl$coded_date,
                                              pheno_yp_nl$REPETITION)
              pheno_c_ynmpp$ntrt <- interaction(pheno_c_ynmpp$PLOT_ID, 
                                                pheno_c_ynmpp$GERMPLASM_ID, 
                                                pheno_c_ynmpp$coded_date,
                                                pheno_c_ynmpp$REPETITION)
              pheno_c_wmffp$ntrt <- interaction(pheno_c_wmffp$PLOT_ID, 
                                                pheno_c_wmffp$GERMPLASM_ID,
                                                pheno_c_wmffp$coded_date,
                                                pheno_c_wmffp$REPETITION)
              names(pheno_yp_nl)[posi] <- 'YNMPL'
              pheno_c_ynmpp$ntrt <- as.character(pheno_c_ynmpp$ntrt)
              pheno_c_wmffp$ntrt <- as.character(pheno_c_wmffp$ntrt)
              pheno_c_wmffp <- pheno_c_wmffp %>% 
                filter(!is.na(ntrt))
              pheno_c_ynmpp <- pheno_c_ynmpp %>% 
                filter(!is.na(ntrt))
              pheno_wp_yp <- pheno_c_wmffp %>% inner_join(pheno_c_ynmpp, by = 'ntrt')
              pheno_yptot <- pheno_wp_yp %>% filter(WP != 0) %>% 
                group_by(PLOT_NUMBER.x) %>% 
                summarise(YPTOT = sum(YP,na.rm = T)) %>% 
                dplyr::select(PLOT_NUMBER.x, YPTOT) %>% 
                right_join(pheno_wp_yp, by = 'PLOT_NUMBER.x') # YPTOT: number of fruits with weight within wach plot
              wt_fr <- as.numeric(pheno_yptot$WP)/as.numeric(pheno_yptot$YPTOT) # weight/fr
              pheno_wp_yp <- pheno_wp_yp[,1:dim(pheno)[2]]
              names(pheno_wp_yp) <- names(pheno)
              names(pheno_wp_yp)[posi] <- 'WTFR'
              pheno_wp_yp$WTFR <- wt_fr
              pheno_yptot_1 <- pheno_yp_nl %>% group_by(PLOT_ID) %>% 
                summarise(YNMPL_TOT = sum(YNMPL, na.rm = T)) %>% 
                right_join(pheno_yp_nl, by = 'PLOT_ID') # YNMPL_TOT: total number of fruit per plot
              
              for (i in 1:length(dates_interval)){
                pheno_wp_yp$coded_date[pheno_wp_yp$OBS_DATE_RECORDED %in% dates_interval[[i]]] <- i
              }
              pheno_wp_yp$ntrt <- interaction(pheno_wp_yp$PLOT_ID, 
                                              pheno_wp_yp$GERMPLASM_ID,
                                              pheno_wp_yp$coded_date,
                                              pheno_wp_yp$REPETITION)
              pheno_wp_yp_1 <- pheno_wp_yp %>% 
                inner_join(pheno_yptot_1, by = 'ntrt')
              WMFPL <- as.numeric(pheno_wp_yp_1$WTFR)*as.numeric(pheno_wp_yp_1$YNMPL_TOT)
              pheno_wp_yp_1 <- pheno_wp_yp_1[,1:dim(pheno)[2]]
              names(pheno_wp_yp_1) <- names(pheno)
              pheno_wp_yp_1$TRAIT_VALUE <- WMFPL
              
              pheno_wp_yp_1$OBSRVTN_REF_CD <- 'WMFPL'
              # pheno_wp_nl$TRAIT_VALUE <- WMFPL
              
              pheno <- rbind(pheno, pheno_wp_yp_1)
              
              # Calculated Trait: AFW_PL = WMFPL/YNMPL
              # average fruit weight per fruit per plot
              pheno_wp_nl <- pheno_wp_yp_1
              pheno_yp_nl <- pheno_yp_nl[,1:dim(pheno)[2]]
              names(pheno_yp_nl) <- names(pheno)
              pheno_yp_nl$ntrt <- interaction(pheno_yp_nl$PLOT_ID, pheno_yp_nl$GERMPLASM_ID)
              pheno_wp_nl$ntrt <- interaction(pheno_wp_nl$PLOT_ID, pheno_wp_nl$GERMPLASM_ID)
              if (length(intersect(pheno_yp_nl$ntrt, pheno_wp_nl$ntrt))>=1){
                ynmpl_sum <- aggregate(TRAIT_VALUE~ntrt, data = pheno_yp_nl, sum)
                pheno_yl_wl <- merge(pheno_wp_nl, ynmpl_sum, by = "ntrt")
                AFW_PL <- pheno_yl_wl$TRAIT_VALUE.x/pheno_yl_wl$TRAIT_VALUE.y
                print(names(pheno_yl_wl))
                print(names(pheno))
                pheno_yl_wl <- pheno_yl_wl[,2:(dim(pheno)[2]+1)]
                names(pheno_yl_wl) <- names(pheno)
                pheno_yl_wl$OBSRVTN_REF_CD <- 'AFW_PL'
                pheno_yl_wl$TRAIT_VALUE <- AFW_PL
                
                pheno <- rbind(pheno, pheno_yl_wl)
              }
            } # end of NHPLT
          }# end of WMFFP
        } # end of case 2
      } # end of NHPLT
    } # prog
  } # crop
  # Calcualte trait for Carrot, 7M, FLBOL
  if (crop == "Carrot"){
    if (prog == '7M'){
      if ('FLBOL' %in% unique(pheno$OBSRVTN_REF_CD)){
        pheno_flbol <- pheno %>% filter(OBSRVTN_REF_CD == 'FLBOL')
        rep_last <- aggregate(REPETITION ~ FIELD_NAME, data = pheno_flbol, max)
        if (nrow(rep_last)>0){
          flbol_last <- pheno_flbol %>% 
            filter(FIELD_NAME == rep_last[1,1] & REPETITION == rep_last[1,2])
          if (nrow(rep_last)>1){
            for (i in 2:nrow(rep_last)){
              temp_last <- pheno_flbol %>% 
                filter(FIELD_NAME == rep_last[i,1] & REPETITION == rep_last[i,2])
              flbol_last <- rbind(flbol_last, temp_last)
            }
          }else{
            flbol_last <- flbol_last
          }
        }
        flbol_last <- flbol_last %>% 
          mutate(OBSRVTN_REF_CD = "FLBOL_LAST")
        pheno <- rbind(pheno, flbol_last)
      }
    }
  }
  
  #   # Calculate trait for Processing Tomato, 9Z, PYLDPA
  if (crop == 'Tomato'){
    if (prog == '9Z'){
      if ('PYLDPA' %in% unique(pheno$OBSRVTN_REF_CD)){
        pheno_pyldpa<- pheno %>% 
          filter(OBSRVTN_REF_CD == 'PYLDPA') %>% 
          mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE) * 0.4047,
                 OBSRVTN_REF_CD = 'PYLDPA_C',
                 UOM = 'tons/acre')
        pheno <- rbind(pheno, pheno_pyldpa)
      }
      if ('YLDPR' %in% unique(pheno$OBSRVTN_REF_CD)){
        pheno_yldpr<- pheno %>% 
          filter(OBSRVTN_REF_CD == 'YLDPR') %>% 
          mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE) * 0.4047,
                 OBSRVTN_REF_CD = 'YLDPR_C',
                 UOM = 'tons/acre')
        pheno <- rbind(pheno, pheno_yldpr)
      }
    }
  }
  pheno
} ## END CALC_TRT FUNCTION


#Interactive functon for repeated measures.
repps= function(pheno){
  vals=list()
  crp=unique(pheno$CROP)
  npheno=subset(pheno, REPETITION > 1)
  trts = unique(npheno$OBSRVTN_REF_CD)
  npheno=subset(pheno,OBSRVTN_REF_CD %in% trts)
  repls=unique(interaction(npheno$OBSRVTN_REF_CD,npheno$REPETITION))
  trt=unlist(lapply(strsplit(as.character(repls),'.',fixed=T),function(x)x[1]))
  rp=unlist(lapply(strsplit(as.character(repls),'.',fixed=T),function(x)x[2]))
  if (dim(npheno)[1]==0) {
    DT=data.table('Crop'=crp,'Trait'=NA,'Repetition'=NA,'TraitXrepetition'='No Repeated Measures')
  }else {
    DT=data.table('Crop'=crp,'Trait'=trt,'Repetition'=rp,'TraitXrepetition'=repls)
  }
  vals[["Data"]]=DT
  vals[["Pheno"]]=pheno
  vals  
}     



Avg_row_head<-function(vals){
  #Delete selected rows after addition
  #vals=vals()
  #rows_to_del=as.numeric(gsub("Row","",input$checked_rows))
  #dat=vals[["Data"]][rows_to_del]
  dat=vals[["Data"]]
  #print(dat)
  #vals$Data=vals$Data[-row_to_del]
  trait=unique(dat$Trait)
  print(trait)
  for (j in trait){
    pheno=vals[["Pheno"]]
    #if (click_value() == 1 & length(rows_to_del) == dim(vals$Data)[1]){
    #  dat <- vals()$Data
    #}else{
    #  dat <- dat
    #}
    
    dat = dat
    #print(dim(pheno))
    dat_j=subset(dat,Trait==j)
    print(dat_j)
    #pheno_n=pheno %>% filter(!OBSRVTN_REF_CD %in% j)
    #print(dim(pheno_n))
    repp=dat_j$Repetition
    print(repp)
    pheno_d=subset(pheno,OBSRVTN_REF_CD == j)
    print(dim(pheno_d))
    pheno_d=subset(pheno_d,REPETITION %in% repp)
    print(dim(pheno_d))
    pheno_d$ID = paste(pheno_d$PEDIGREE_NAME, pheno_d$REP_NUMBER, pheno_d$OBSRVTN_REF_CD, pheno_d$GROWSEASON, 
                       pheno_d$TEST_SET_NAME, pheno_d$FIELD_NAME)
    pheno_d$TRAIT_VALUE = as.numeric(pheno_d$TRAIT_VALUE)
    a = duplicated(pheno_d$ID)
    ndat = pheno_d[!a,]
    pheno_d$OBSRVTN_REF_CD=as.character(pheno_d$OBSRVTN_REF_CD)
    Mean = tapply(pheno_d$TRAIT_VALUE, list(pheno_d$ID,pheno_d$OBSRVTN_REF_CD), mean,na.rm=T)
    #for (i in 1:ncol(Mean)){
    print(head(Mean))
    Mean_val = sapply(1:length(Mean), function(k) Mean[[k]])
    L = match(ndat$ID, rownames(Mean))
    ndat$TRAIT_VALUE = Mean_val[L]
    ndat$REPETITION = paste0("RepAvg_",paste(min(as.numeric(repp)),max(as.numeric(repp)),sep='-'))
    ndat$OBSRVTN_REF_CD=paste0(j,"Avg_",paste(min(as.numeric(repp)),max(as.numeric(repp)),sep='-'))
    ndat=subset(ndat,select=-ID)
    npheno=rbind(pheno,ndat)
    print(unique(npheno$OBSRVTN_REF_CD))
    vals[["Pheno"]]=npheno
    #addedtraits <- paste(newtraits(),unique(ndat$OBSRVTN_REF_CD), sep = "\n")
    #newtraits(addedtraits)
  }
  return(vals)
}

#####
## Sum
Sum_rep<-function(vals){
  dat=vals[["Data"]]
  #rows_to_del=as.numeric(gsub("Row","",input$checked_rows))
  #dat=vals$Data[rows_to_del]
  #vals$Data=vals$Data[-row_to_del]
  trait=unique(dat$Trait)
  for (j in trait){

    pheno=vals[["Pheno"]]
    
    #  pheno=vals$Pheno
    #  if (click_value() == 1 & length(rows_to_del) == dim(vals$Data)[1]){
    #    dat <- vals()$Data
    #  }else{
    #    dat <- dat
    #  }
    #pheno_n=pheno %>% filter(!OBSRVTN_REF_CD %in% j)
    dat_j=subset(dat,Trait==j)
    repp=dat_j$Repetition
    pheno_d=subset(pheno,OBSRVTN_REF_CD %in% j)
    pheno_d=subset(pheno_d,REPETITION %in% repp)
    #### Temporary fix for Cucumber calculated traits NFSA_C and WFSA_C #########

      ndat <- pheno_d
      ndat_1 <- ndat %>% 
        # select( TRAIT_VALUE, ASORT1, FIELD_NAME, REPETITION, PLOT_BID, ABS_C, ABS_R) %>% 
        dplyr::group_by(CROSS_NAME,FIELD_NAME,PLOT_BID,ABS_C,ABS_R,REPETITION, GROWSEASON, OBSRVTN_REF_CD) %>% 
        dplyr::summarize(mean_for_sq = mean(TRAIT_VALUE)) %>% 
        dplyr::ungroup() %>% 
        dplyr::group_by(CROSS_NAME,FIELD_NAME,PLOT_BID,ABS_C,ABS_R, GROWSEASON, OBSRVTN_REF_CD) %>% 
        dplyr::summarize(TRAIT_VALUE = sum(mean_for_sq)) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(OBSRVTN_REF_CD=paste0(j,"Sum_",paste(min(as.numeric(repp)),max(as.numeric(repp)),sep='-')),
               REPETITION = paste0("RepSum_",paste(min(as.numeric(repp)),max(as.numeric(repp)),sep='-')))
      
      print("Cucumber repeated measure")
      print(head(ndat_1))
      print('dimension of ndat_1 and ndat')
      print(dim(ndat_1))
      print(dim(ndat))
      
      ndat <- ndat %>%
        dplyr::select(-OBSRVTN_REF_CD) %>%
        dplyr::select(-REPETITION) %>%
        dplyr::select(-TRAIT_VALUE) %>%
        full_join(ndat_1, by = c('CROSS_NAME','FIELD_NAME','PLOT_BID','ABS_C','ABS_R', 'GROWSEASON'))
      
      ndat <- ndat[!duplicated(ndat[, c('CROSS_NAME', 'PLOT_BID', 'GROWSEASON', 'FIELD_NAME')]),]
      
      print('ndat after aggregation')
      print(head(ndat))
      
      print('number of cn')
      print(length(unique(ndat$CROSS_NAME)))
      
      ndat <- ndat[,colnames(pheno_d)%in% colnames(ndat) ]

    npheno=rbind(pheno,ndat)
    print(unique(npheno$OBSRVTN_REF_CD))
    vals[["Pheno"]]=npheno
    #addedtraits <- paste(newtraits(),unique(ndat$OBSRVTN_REF_CD), sep = "\n")
    #newtraits(addedtraits)
  }
  
  return(vals)
  
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

