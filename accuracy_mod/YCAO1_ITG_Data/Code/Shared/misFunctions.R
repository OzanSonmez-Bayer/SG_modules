CompareReliabilityPlot <- function(file1, file2, trait_name){
  
  # file1 and file2 are strings: names of files for BLUP with A matrix and Deeper A matrix respectively
  # trait: a stirng, trait name
  
  old_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file1)))
  old_con <- textConnection(old_data_file)
  old_data <- read.csv(old_con)
  close(old_con)
  new_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file2)))
  new_con <- textConnection(new_data_file)
  new_data <- read.csv(new_con)
  close(new_con)
  
  overlapping_ped <- intersect(as.character(old_data$PEDIGREE_NAME), as.character(new_data$PEDIGREE_NAME))
  reliability_dat <- rbind(new_data, old_data)
  reliability_dat <- reliability_dat %>% 
    mutate(Amatrix = as.factor(rep(c('5 Gen', "Parent Only"), c(nrow(new_data), nrow(old_data))))) %>% 
    filter(trait == trait_name & PEDIGREE_NAME %in% overlapping_ped)
  
  
  p <- ggplot(reliability_dat, aes(x = Amatrix, y = reliability, fill = Amatrix)) + 
    geom_violin() + 
    stat_summary(fun.y = mean, geom = 'point', color = "blue", size = 5) +
    labs(title = paste("Reliability: ", trait_name)) + 
    xlab("A Matrix")
  p
}


CompareBLUP <- function(file1, file2, trait_name){
  
  # file1 and file2 are strings: names of files for BLUP with A matrix and Deeper A matrix respectively
  # trait: a stirng, trait name
  
  old_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file1)))
  old_con <- textConnection(old_data_file)
  old_data <- read.csv(old_con)
  close(old_con)
  new_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file2)))
  new_con <- textConnection(new_data_file)
  new_data <- read.csv(new_con)
  close(new_con)
  
  overlapping_ped <- intersect(as.character(old_data$PEDIGREE_NAME), as.character(new_data$PEDIGREE_NAME))
  new_data_1 <- new_data %>% 
    filter(as.character(PEDIGREE_NAME) %in% overlapping_ped & trait == trait_name) %>% 
    dplyr::select(PEDIGREE_NAME, predicted.value) %>% 
    dplyr::rename(five_gen_blup = predicted.value)
  old_data_1 <- old_data %>% 
    filter(as.character(PEDIGREE_NAME) %in% overlapping_ped & trait == trait_name) %>% 
    dplyr::select(PEDIGREE_NAME, predicted.value) %>% 
    dplyr::rename(two_gen_blup = predicted.value)
  
  BLUP_dat <- new_data_1 %>% 
    join(old_data_1, by = 'PEDIGREE_NAME')

  
  # p <- ggplot(BLUP_dat, aes(x = two_gen_blup, y = five_gen_blup)) + 
  #   geom_point(alpha = 0.5) + 
  #   annotate(geom = 'text',
  #            label = paste('blup cor = ', round(cor(BLUP_dat$two_gen_blup, BLUP_dat$five_gen_blup, use = "complete.obs", 
  #                                                   method = 'pearson'),2)),
  #            x = min(BLUP_dat$two_gen_blup, na.rm = T) + 0.8,
  #            y = max(BLUP_dat$five_gen_blup, na.rm = T) - 0.3) + 
  #   annotate(geom = 'text',
  #            label = paste('rank cor = ', round(cor(BLUP_dat$two_gen_blup, BLUP_dat$five_gen_blup, use = "complete.obs",
  #                                                   method = 'spearman'),2)),
  #            x = min(BLUP_dat$two_gen_blup, na.rm = T) + 0.8,
  #            y = max(BLUP_dat$five_gen_blup, na.rm = T) - 0.9) + 
  #   labs(title = paste('BLUP: ', trait_name), x = 'Parent Only BLUP', y = 'Five Gen Acestors BLUP') + 
  #   theme_minimal()
  p <- ggplot(BLUP_dat, aes(x = two_gen_blup, y = five_gen_blup)) + 
    geom_point(alpha = 0.5, color = 'blue') + 
    labs(title = paste('BLUP: ', trait_name), x = 'Parent Only BLUP', y = 'Five Gen Acestors BLUP') + 
    theme_minimal() + 
    theme(text = element_text(size = 15))
  p
}


CompareCorrelation <- function(file1, file2, trt_list){
  
  # file1 and file2 are strings: names of files for BLUP with A matrix and Deeper A matrix respectively
  # trt_list: a list of trait of interests
  rank_df <- data.frame()
  for (i in 1:length(trt_list)){
    old_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file1)))
    old_con <- textConnection(old_data_file)
    old_data <- read.csv(old_con)
    close(old_con)
    new_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file2)))
    new_con <- textConnection(new_data_file)
    new_data <- read.csv(new_con)
    close(new_con)
    
    overlapping_ped <- intersect(as.character(old_data$PEDIGREE_NAME), as.character(new_data$PEDIGREE_NAME))
    new_data_1 <- new_data %>% 
      filter(as.character(PEDIGREE_NAME) %in% overlapping_ped & trait == trt_list[i]) %>% 
      dplyr::select(PEDIGREE_NAME, predicted.value) %>% 
      dplyr::rename(five_gen_blup = predicted.value)
    old_data_1 <- old_data %>% 
      filter(as.character(PEDIGREE_NAME) %in% overlapping_ped & trait == trt_list[i]) %>% 
      dplyr::select(PEDIGREE_NAME, predicted.value) %>% 
      dplyr::rename(two_gen_blup = predicted.value)
    
    BLUP_dat <- new_data_1 %>% 
      join(old_data_1, by = 'PEDIGREE_NAME')
    
    temp_cor <- data.frame('Trait' = trt_list[i],
                            'BLUP Correlation' = round(cor(BLUP_dat$two_gen_blup, BLUP_dat$five_gen_blup,
                                                    use = 'complete.obs',
                                                    method = 'pearson'),2),
                           'Rank Correlation' = round(cor(BLUP_dat$two_gen_blup, BLUP_dat$five_gen_blup,
                                                          use = 'complete.obs',
                                                          method = 'spearman'),2))
    if (i == 1){
      rank_df <- temp_cor 
    }else{
      rank_df <- rbind(rank_df, temp_cor)
    }
  }
  return(rank_df)
}

visualCorr <- function(dataset, title_name){
  # dataset: a data frame with correlation
  # title_name: a string
  p <- ggplot(dataset, aes(x = Trait, y = BLUP.Correlation, color = mdl)) + 
    geom_point(shape = 1, size = 4, stroke = 2) + 
    theme_minimal() + 
    theme(text = element_text(size = 13),
          axis.text.x = element_text(size = 13, angle = 90)) +
    labs(title = title_name)
  p
}

heritabilityVis <- function(file1, file2, title_name){
  # file1: h2 from 2-gen A matrix
  # file2: h2 from 5-gen A matrix
  # title_name: plot title
  old_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file1)))
  old_con <- textConnection(old_data_file)
  old_data <- read.csv(old_con)
  close(old_con)
  new_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file2)))
  new_con <- textConnection(new_data_file)
  new_data <- read.csv(new_con)
  close(new_con)
  
  old_data <- old_data %>% 
    dplyr::rename(Parent_only = h2) %>% 
    dplyr::select(trait, Parent_only)
  new_data <- new_data %>% 
    dplyr::rename(Five_Gen_Ancestor = h2) %>% 
    dplyr::select(trait,Five_Gen_Ancestor)
  dataset <- old_data %>% 
    join(new_data, by = 'trait')
  set.seed(2019)
  p <- ggplot(dataset, aes(x = Parent_only, y = Five_Gen_Ancestor, label = trait)) + 
    geom_point(color = 'blue') + 
    geom_text_repel(aes(label = trait), size = 3.5) + 
    geom_abline(intercept = 0, slope = 1) + 
    theme_minimal() + 
    labs(title = title_name)
  p
}

# Check lines without direct observed hybrid data

linesCheck <- function(file1){
  # file1: a data frame with PGCA results
  gca_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file1)))
  con_file <- textConnection(gca_file)
  gca_dat <- read.csv(con_file)
  close(con_file)
  unique_line <- length(unique(gca_dat$PEDIGREE_NAME))
  unique_new_line <- length(unique(gca_dat$PEDIGREE_NAME[which(is.na(gca_dat$n))]))
  return(list('total_line' = unique_line, 'new_line' = unique_new_line, 'percent' = round(unique_new_line/unique_line,2) * 100))
}

# Function CompareReliabilityPlot_other is to check the reliability of other relatives
CompareReliabilityPlot_other <- function(file1, file2, trait_name){
  
  # file1 and file2 are strings: names of files for BLUP with A matrix and Deeper A matrix respectively
  # trait: a stirng, trait name
  
  old_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file1)))
  old_con <- textConnection(old_data_file)
  old_data <- read.csv(old_con)
  close(old_con)
  new_data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', file2)))
  new_con <- textConnection(new_data_file)
  new_data <- read.csv(new_con)
  close(new_con)
  
  overlapping_ped <- intersect(as.character(old_data$PEDIGREE_NAME), as.character(new_data$PEDIGREE_NAME))
  reliability_dat <- new_data %>% 
    filter(trait == trait_name) %>% 
    filter(!as.character(PEDIGREE_NAME) %in% overlapping_ped) %>% 
    mutate(mdl = '5 Gen')
  
  
  p <- ggplot(reliability_dat, aes(x = mdl, y = reliability, fill = mdl)) + 
    geom_violin() + 
    stat_summary(fun.y = mean, geom = 'point', color = "blue", size = 5) +
    labs(title = paste("Reliability of Other Relatives: ", trait_name))
  p
}

vizBLUP <- function(blup_df, blup_type, trait_name){
  # blup_df: a string, file name which has results from deep pedigree
  # blup_type: a string, P_SCA or P_GCA
  # trait_name: a string, trait name
  
  data_file <- rawToChar(get_object(paste0('s3://','genome-analytics-perm-space/', blup_df)))
  con <- textConnection(data_file)
  blup_df <- read.csv(con)
  close(con)
  
  blup_sub <- blup_df %>% 
    filter(trait == trait_name & mdl == blup_type) %>% 
    mutate(is_direct = if_else(is.na(n), 0, 1)) %>% 
    mutate(is_direct = as.factor(is_direct)) %>% 
    select(predicted.value, is_direct)
  
  p <- ggplot(blup_sub, aes(y = predicted.value, x = is_direct)) +
    geom_point(aes(color = is_direct))
  p
  
}

cal_trt <- function(pheno){
  pheno$TRAIT_VALUE=as.numeric(pheno$TRAIT_VALUE)
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
    NFSA_C=ft_mk/n_pt*n_pa
    
    WFSA_C=nt_wt/n_pt*n_pa
    pheno_comb=pheno_comb[,2:(dim(pheno)[2]+1)]
    names(pheno_comb)=names(pheno)
    
    nt_dat <- pheno_c_nt %>% 
      group_by(PLOT_BID) %>% 
      dplyr::summarise(sum_nt = sum(NT, rm.na = T))
    
    
    fk_dat <- pheno_c_fk %>% 
      group_by(PLOT_BID) %>% 
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
  return(pheno)
}


outL = function(X){
  if (length(X)==0){out=NULL}else{
    pval = grubbs.test(X)$p.value
    if (pval>0.01 | is.na(pval)){ out=X }else{
      out = rm.outlier(X, fill=TRUE,median=TRUE)
    }
  }
  list(out = out, pval = pval)
}

rm.out = function(pheno){
  trait = unique(pheno$OBSRVTN_REF_CD)
  test = unique(pheno$TEST_SET_NAME)
  pheno$TRAIT_VALUE=as.numeric(as.character(pheno$TRAIT_VALUE))
  #data.new = pheno
  Tab = data.frame(Row=c(),Cross = c(), SetName = c(),PlotBid=c(), Trait =  c(), Season =  c(),
                   Pvalue =  c(), TraitVal = c(), mean =  c(), 
                   median =  c(), SD =  c())
  
  for (i in 1:length(trait)){
    for (j in 1: length(test)){
      dat = subset(pheno, TEST_SET_NAME==test[j] & OBSRVTN_REF_CD==trait[i])
      if (dim(dat)[1]!=0){
        
        tesT = outL(as.numeric(dat$TRAIT_VALUE))
        c = tesT$out
        p = tesT$pval

        if (sum(is.na(c))==length(c)) c=as.numeric(dat$TRAIT_VALUE)
        out.w = which(c!=dat$TRAIT_VALUE)
        m = length(out.w)
        data.out = dat[out.w, ]
        Mean = round(mean(dat$TRAIT_VALUE),2)
        Median = round(median(dat$TRAIT_VALUE),2)
        SD = round(sd(dat$TRAIT_VALUE),2)
        # Sample = paste(dat$TRAIT_VALUE, collapse = ",")
        # Series=paste(hist(dat$TRAIT_VALUE,plot=F)$counts,collapse=",")
        V = data.frame(Row=rownames(data.out),
                       Cross=data.out$CROSS_NAME, 
                       SetName = rep(test[j], m), 
                       PlotBid=data.out$PLOT_BID,
                       Obs_Date_Recorded = data.out$OBS_DATE_RECORDED,
                       Obs_Date_Modified =  data.out$OBS_DATE_MODIFIED,
                       Trait = rep(trait[i], m),
                       Season = data.out$GROWSEASON,
                       Pvalue = rep(round(p,5), m), 
                       TraitVal = round(data.out$TRAIT_VALUE,2), 
                       mean = rep(Mean, m), 
                       median = rep(Median, m), 
                       SD = rep(SD,m))
        Tab = rbind(Tab, V)
      }
    }
  }
  Tab
}

pullData <- function(season, setlist, trtlist, crop, selectedProg){
  # season: a vector of strings. i.e., c('2019:01', '2019:02')
  # setlist: a vector of string. set names
  # trtlist: a vector of stirng. traits of interests
  # crop: a string.
  # selectedProg: a string. Program name
  
  source('Credential/Credential_h2h.R')
  path <- '/shared/Veg_Phenotypes/H2H_Data/'
  
  data_df <- data.frame()
  for (SelectedSeason in season){
    print(SelectedSeason)
    SelectedSeason <- gsub(':', '/', SelectedSeason)
    infile <- paste0(path,crop,"/",SelectedSeason,"/",selectedProg,".RData")
    
    if(head_object(object = infile, bucket = "genome-analytics-perm-space")[[1]]){
      print("OBJECT FOUND")
      s3load(object = infile, bucket = "genome-analytics-perm-space")
      Phenos <- Phenos %>% 
        filter(TEST_SET_NAME %in% setlist & OBSRVTN_REF_CD %in% trtlist)
      data_df <- plyr::rbind.fill(data_df, Phenos)
    } else {
      print("OBJECT NOT FOUND")
    }
  }
  data_df$CROSS_NAME[which(data_df$CROSS_NAME == "")] <- data_df[which(data_df$CROSS_NAME == ''), 'PEDIGREE_NAME']
  data_df$TREP <- paste(data_df$REP_NUMBER,data_df$TRACK_ID, sep = '_')
  
  return(data_df)
}

# Check Unit of Measure
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
          }else{
            #Generic function
            dat$TRAIT_VALUE=as.numeric(dat$TRAIT_VALUE)
            mn1=mean(as.numeric(dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE),na.rm=T)
            mn2=mean(as.numeric(dat[dat$UOM==as.character(names(tab)[2]),]$TRAIT_VALUE),na.rm=T)
            rtio=mn1/mn2
            if (round(rtio)>2) {
              conv=round(dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE/rtio,3)
            }else{
              conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
            }
          } # END IF UNITS NOT FOUND
          
        }else{
          conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
        } # END IF RATION OF UNITS > 3
      }else{
        conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
      } # end if tmp 1 or tmp 2 is NA
    }else{
      conv=dat[dat$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE
    } # END IF THERE'S ONLY ONE UNIT
    pheno[pheno$OBSRVTN_REF_CD==trait[i] & pheno$UOM==as.character(names(tab)[1]),]$TRAIT_VALUE=conv
  }
  pheno
}


hrFun <- function(dat){
  # dat: data frame. SCA results return by function pblup_fun
  # reutrn: data frame. trait, heritability, average reliability per trait, and the number of observations
  sum_tab <- dat %>% 
    group_by(trait) %>% 
    dplyr::summarise(nline = length(unique(PEDIGREE_NAME)),
              h2 = unique(h2),
              avg_r2 = round(mean(reliability),3))
  return(sum_tab)
}


AWS.Authorization <- function (username) {
  # username: ycao1 or sdi
  if (username == 'ycao1'){
    VaultPathToClientSecret = "/secret/veg-pipeline-solutions/prod/aws/ycao1"
    VAULT_PATH <- VaultPathToClientSecret
    role_id="4e9c1e37-a380-2b6a-92d7-0b6add761f00"
    secret_id=" beb985d2-70ed-02c8-e24d-56899dcd972b"
    secret_id_accessor="6bc5584f-70ec-ec8d-d2c0-6cbfb25a313c" 
  }
  submit_body_list <- toJSON(list(role_id=role_id, secret_id=secret_id), auto_unbox = T)
  vault_client_secret <- function(secret_path, submit_bodyList) {
    Base_URL <- "https://vault.agro.services/v1"
    VAULT_BASE_URL <- paste(Base_URL,"/auth/approle/login", sep = "")
    url <- VAULT_BASE_URL
    response <- httr::POST(url,
                           body = submit_bodyList,
                           encode = "json")
    #print(response)
    VAULT_TOKEN <- httr::content(response, as = "parsed", encoding = "UTF-8")$auth$client_token
    #print(VAULT_TOKEN)
    VAULT_PATH_URL <- paste(Base_URL, secret_path, sep = "")
    vault_list <- httr::GET(VAULT_PATH_URL,
                            httr::add_headers("X-Vault-Token" = VAULT_TOKEN)
    )
    vault_secret <- httr::content(vault_list, as = "parsed", encoding = "UTF-8")
    keyID <-vault_secret$data$AWS_ACCESS_KEY_ID
    secret <-vault_secret$data$AWS_SECRET_ACCESS_KEY
    Sys.setenv(AWS_ACCESS_KEY_ID = keyID)
    Sys.setenv(AWS_SECRET_ACCESS_KEY = secret)
  }
  PING_URL <- "https://amp.monsanto.com/as/token.oauth2"
  vault_client_secret(VAULT_PATH, submit_body_list)
}


repeatedMeasure <- function(dat, trait, season, mdl){
  # dat: data frame
  # trait: a string
  # mdl: a string, aggregation method, sum or average
  
  dat_trt <- dat %>% 
    filter(OBSRVTN_REF_CD == trait & GROWSEASON == season) %>% 
    mutate(ID  = paste(PEDIGREE_NAME, REP_NUMBER, TEST_SET_NAME, FIELD_NAME, ENTRY_NUM),
           TRAIT_VALUE = as.numeric(TRAIT_VALUE))
  repp <- as.numeric(unique(dat_trt$REPETITION))
  
  ndat <- dat_trt
  # ndat <- dat_trt[!duplicated(dat_trt$ID), ]
  if (mdl %in% c('mean', 'sum')){
    if (mdl == 'mean'){
      calc_val <- tapply(dat_trt$TRAIT_VALUE, list(dat_trt$ID), mean, na.rm = T)
      calc_val_1 <- sapply(1:length(calc_val), function(k){calc_val[[k]]})
    }else if (mdl == 'sum'){
      calc_val <- tapply(dat_trt$TRAIT_VALUE, list(dat_trt$ID), sum, na.rm = T)
      calc_val_1 <- sapply(1:length(calc_val), function(k){calc_val[[k]]})
    }
    
    L <- match(ndat$ID, rownames(calc_val))
    ndat$TRAIT_VALUE <- calc_val_1[L]
    ndat$REPETITION <- 1
    ndat$OBSRVTN_REF_CD <- paste0(trait, '_', mdl, "_", paste(min(repp), max(repp), sep = "-"))
    
    ndat <- subset(ndat, select = -ID) 
  }else{
    ndat <- dat_trt
  }
  return(ndat)
}

writing_csv <- function(dat, filename){
  write.csv(dat, file = filename, row.names = F)
}

calc_stage_count <- function(dat){
  dat_female <- dat %>% 
    dplyr::select(PEDIGREE_NAME, P1, EXPER_STAGE_REF_ID) %>% 
    dplyr::rename(pedigree_name = PEDIGREE_NAME, 
                  parent_pedigree_name = P1, 
                  stage = EXPER_STAGE_REF_ID)
  
  dat_male <- dat %>% 
    dplyr::select(PEDIGREE_NAME, P2, EXPER_STAGE_REF_ID) %>% 
    dplyr::rename(pedigree_name = PEDIGREE_NAME, 
                  parent_pedigree_name = P2, 
                  stage = EXPER_STAGE_REF_ID)
  
  dat_concat <- rbind.data.frame(dat_female, dat_male)
  
  stage_count_df = dat_concat %>% 
    dplyr::filter(parent_pedigree_name != '') %>% 
    dplyr::group_by(parent_pedigree_name, stage) %>% 
    dplyr::summarize(stage_count = length(unique(pedigree_name)))
  
  return(stage_count_df)
}
