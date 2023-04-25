# Geno functions

filterfun <- function(geno){
  # Remove individuals with more than a certain % missing data
  individualmissing <- apply(geno, 1, function(x){
    return(length(which(is.na(x)))/ncol(geno))
  })
  
  # Remove markers with certain % missing data
  markermissing <- apply(geno, 2, function(x){
    return(length(which(is.na(x)))/nrow(geno))
  })
  return(list(individualmissing, markermissing))
}

# Function getAlleles and gen2Additive from genoPull
# Convert ACGT to 0, 1, 2
getAlleles=function(gen){
  allele=lapply(apply(gen,2,unique),function(x)sort(na.omit(unique(unlist(strsplit(x,''))))))
  
  #$NCSAT009225482
  #[1] ";" "G" "T"	#If uncollapsed pedigrees will cause ; to still be in there...
  #$NCSAT009220828
  #[1] "G" 			#For single allele
  
  #gen$allele=allele
  #gen=gen[nchar(gen$allele)!=1,]		#Don't think I want to filter anything just return gen with allele column 
  #gen=gen[is.na(gen$allele) != 1,]
  allele
}

gen2Additive <- function( gen ) {
  allele=getAlleles(gen)
  #Check for markers with >2 alleles: exclude them
  allen=lapply(allele,function(x)length(x))
  #Check for markers where alleles are not (ACTG): exclude them
  wrg_al=lapply(lapply(allele,function(x)!x %in% c('A','C','G','T')),sum)
  #excl.markers=names(which(allen !=2))	#note markers with <> 2 alleles
  #ngen=matrix(5,nrow=nrow(gen),ncol=ncol(gen))
  
  #Only keep markers with 1 or 2 allele	#Modified to keep 3/20 fixed markers. Add allen==2 for polymorphic markers.
  gd_mrkr=names(which(allen >=1 & allen <=2))			
  #Only keep markers with alleles (ACTG)
  gd_mrkr=gd_mrkr[! gd_mrkr %in% names(which(wrg_al >= 1))]	#Remove all Indels ***,Ins,Del etc. Only accepts bp.
  old=ncol(gen)
  print(c('Trouble Markers',allele[wrg_al >=1]))
  gen <- gen[,is.element(colnames(gen),gd_mrkr)]
  allele = allele[names(allele) %in% gd_mrkr]
  new=ncol(gen)
  print(c('#Markers removed for having *, NULL | > 2 alleles:',old-new))
  #
  for(i in 1:ncol(gen)) {
    alleles=c(paste(allele[[i]][1],allele[[i]][1],sep=''),paste(allele[[i]][1],allele[[i]][2],sep=''),paste(allele[[i]][2],allele[[i]][2],sep=''))
    if (alleles[1] %in% gen[,i]){
      gen[,i]=sub(alleles[1],as.numeric(2),gen[,i])}
    if (alleles[2] %in% gen[,i]){
      gen[,i]=sub(alleles[2],as.numeric(1),gen[,i])}
    if (alleles[3] %in% gen[,i]){
      gen[,i]=sub(alleles[3],as.numeric(0),gen[,i])}
  }
  #ngen[is.na(gen)]=5
  sname=colnames(gen)
  gen=apply(gen,2,as.numeric)	#Makes sure numeric
  sname=sub('Q-','',sname)	#Makes sure no 'technology codes included. May need to expand.
  colnames(gen)=sname
  
  gen
}

genClean <- function(gen){
  allele=getAlleles(gen)
  #Check for markers with >2 alleles: exclude them
  allen=lapply(allele,function(x)length(x))
  #Check for markers where alleles are not (ACTG): exclude them
  wrg_al=lapply(lapply(allele,function(x)!x %in% c('A','C','G','T')),sum)
  #excl.markers=names(which(allen !=2))	#note markers with <> 2 alleles
  #ngen=matrix(5,nrow=nrow(gen),ncol=ncol(gen))
  
  #Only keep markers with 1 or 2 allele	#Modified to keep 3/20 fixed markers. Add allen==2 for polymorphic markers.
  gd_mrkr=names(which(allen >=1 & allen <=2))			
  #Only keep markers with alleles (ACTG)
  gd_mrkr=gd_mrkr[! gd_mrkr %in% names(which(wrg_al >= 1))]	#Remove all Indels ***,Ins,Del etc. Only accepts bp.
  old=ncol(gen)
  print(c('Trouble Markers',allele[wrg_al >=1]))
  gen <- gen[,is.element(colnames(gen),gd_mrkr)]
  allele = allele[names(allele) %in% gd_mrkr]
  new=ncol(gen)
  print(c('#Markers removed for having *, NULL | > 2 alleles:',old-new))
  #ngen[is.na(gen)]=5
  sname=colnames(gen)
  gen=apply(gen,2,as.character)	#Makes sure numeric
  sname=sub('Q-','',sname)	#Makes sure no 'technology codes included. May need to expand.
  colnames(gen)=sname
  
  gen
}

genElement <- function(x){
  for (i in 1: ncol(x)){
    for (j in 1:nrow(x)){
      x[j,i] <- paste(strsplit(x[j,i], "|")[[1]][strsplit(x[j,i], "|")[[1]] != '|'], collapse = '')
      if (length(strsplit(x[j, i],split = '')[[1]]) > 2){
        x[j, i] <- NA
      }
    }
  }
  x
}

get_germID <- function(pedigree){ # one pedigree may have more than one germplasm id (!)
  id = unique(CropIDs$M.GERMPLASM.X_ID[CropIDs$M.GERMPLASM.PEDIGREE == pedigree])
  if (length(id) == 0){
    return(data.frame(pedigree = pedigree,
                      germplasm_id = NA,
                      n = 0,
                      stringsAsFactors = F))
    print(paste("pedigree not found:", pedigree))
  } else{
    return(data.frame(pedigree = pedigree,
                      germplasm_id = id,
                      n = 1:length(id),
                      stringsAsFactors = F))
  }
}

KinshipMatrix <- function(Gmatrix){
  # Gmatrix: marker based relateness matrix (weighted or unweighted)
  # res: a data frame, G inverse in Asreml readable format
  rMinv <- ginv(Gmatrix)
  
  res <- data.frame(Row = rep(1:nrow(rMinv), nrow(rMinv)),
                    Column = rep(1:nrow(rMinv), each = nrow(rMinv)),
                    coeff = as.numeric(rMinv), 
                    lower = as.logical(lower.tri(rMinv, diag = TRUE)))
  rm(rMinv)
  # only use lower triangle
  
  res <- res[res$lower == TRUE, c('Row', 'Column', 'coeff')]
  res <- res[order(res$Row, res$Column),]
  res <- res[res$coeff != 0,]
  
  return(res)
  
}

unweightedGFun <- function(geno_file, pheno_file, map_file, plot_draw = TRUE){
  gp <- create.gpData(geno = geno_file, pheno = pheno_file, pedigree = NULL, map = map_file)
  summary(gp)
  # plot(gp$map)
  
  ##### Impute Missing Data
  gp_num <- codeGeno(gp, label.heter = 'alleleCoding',
                     maf = 0.05,
                     nmiss = 0.25,
                     impute = TRUE,
                     impute.type = 'random',
                     verbose = TRUE)
  
  ##### Method 1: VanRaden, G_Base
  
  realized <- kin(gp_num, ret = 'realized')
  summary(realized)
  
  if(summary(realized)$dim['nrow'] > summary(realized)$rank){
  # Remove singularity in the matrix
  
    realizedPD <- nearPD(realized, keepDiag = T)
    G_base <- matrix(realizedPD[[1]]@x, nrow = realizedPD[[1]]@Dim[1])
    G_base <- G_base + diag(0.01, nrow(G_base))
    attr(G_base, 'dimnames') <- realizedPD[[1]]@Dimnames
    class(G_base) <- 'relationshipMatrix'
    str(G_base)
    summary(G_base)
  }else{
    G_base <- realized
  }
  # dimension                     125 x 125 
  # rank                          125 
  # range of off-diagonal values  -0.5662422 -- 2.028777 
  # mean off-diagonal values      -0.01548537 
  # range of diagonal values      0.9796846 -- 2.800697 
  # mean diagonal values          1.930186 
  if (plot_draw){
    plot(G_base)
  }
  return(list(G = G_base,
              gp_obj = gp_num))
}

weightedG_1 <- function(markerM){
  # weightedG_1 function:
  # W = 2*p*(1-p)
  # markers contribute to genomic relatedness proportional to the reciprocal of their expected variance
  # markerM: marker matrix with 0, 1, 2
  # return: G_base, a weighted matrix
  
  allele_p <- apply(markerM, 2, mean)/2
  allele_p <- round(allele_p, 3)
  allele_df <- data.frame(marker = names(allele_p), p = unname(allele_p))
  
  # Center Marker matrix using the means (2 times the allele frequencies)
  allele_mat <- matrix(rep(allele_p*2, nrow(markerM)), ncol = ncol(markerM), nrow = nrow(markerM), byrow = TRUE)
  rownames(allele_mat) <- rownames(markerM)
  colnames(allele_mat) <- colnames(markerM)
  
  # Create Z matrix, centered Marker matrix
  
  Z <- markerM - allele_mat
  
  
  # Generic G matrix  calculation; G = Z*D*Z_transpose. D is the scaling factor, SNP weights
  
  
  D_base <- diag(1 / ncol(markerM) * (2 * allele_p * (1 - allele_p)))
  colnames(D_base) <- colnames(markerM)
  rownames(D_base) <- colnames(markerM)
  
  G_base <- Z %*% D_base %*% t(Z)
  
  return(G_base)
  
}


weightedG_BL <- function(phenodata, markerM, r2, nIter, burnIn, outPath){
  # Bayesian LASSO
  # Parameters:
               # phenodata: a vector of phenotypic value
               # markerM: a matrix with values 0, 1, 2
               # r2: numeric, phenotypic variance attributed to residual variance. 1- heritability
               # nIter: integer, number of iterations
               # burnIn: integer, number of burn in from Gibbs Sampler
               # outPath: string, path and file name to save results
  
  Mrk_sub <- markerM - 1
  mode.sigE <- r2*var(phenodata)
  dfe <- 3
  Se <- mode.sigE*(dfe + 2)
  
  # lambda
  mode.sigL <- (1 - r2)*var(phenodata)
  rate <- 2*Se/mode.sigL*sum(colMeans(Mrk_sub)^2)
  shape <- 1.1
  
  # Set priors
  prior <- list(varE = list(S0 = Se, df0 = dfe),
                lambda = list(type = 'random',
                              shape = shape, 
                              rate = rate))
 
  
  ETA <- list(MRK = list(X = Mrk_sub, model = 'BL', prior))
  fmR_l <- BGLR(y = phenodata, 
                ETA = ETA, 
                nIter = nIter, 
                burnIn = burnIn, 
                saveAt = outPath)
  
  mu <- fmR_l$ETA$MRK$b
  
  allele_p <- apply(markerM, 2, mean)/2
  allele_p <- round(allele_p, 3)
  allele_df <- data.frame(marker = names(allele_p), p = unname(allele_p))
  
  # Center Marker matrix using the means (2 times the allele frequencies)
  allele_mat <- matrix(rep(allele_p*2, nrow(markerM)), ncol = ncol(markerM), nrow = nrow(markerM), byrow = TRUE)
  rownames(allele_mat) <- rownames(markerM)
  colnames(allele_mat) <- colnames(markerM)
  
  # Create Z matrix, centered Marker matrix
  
  Z <- markerM - allele_mat
  
  D_BL <- diag(1 / ncol(markerM) * (2 * allele_p * (1 - allele_p) * mu * mu))
  colnames(D_BL) <- colnames(Mrk_sub)
  rownames(D_BL) <- colnames(Mrk_sub)
  
  G_BL <- Z %*% D_BL %*% t(Z)
  return(G_BL)
}

get_varcomps <- function(c){
  vcs <- summary(c)$varcomp$component
  names(vcs) <- c('Vg', 'Vresid')
  return (vcs)
}

get_pred <- function(r, pred.list){
  ID <- as.character(r['pedigree'])
  ID2 <- paste0('giv(pedigree, var = T)_', ID)
  trait <- paste0("Y", r['fold'])
  return(pred.list[[trait]][ID2, ])
}

missingInFold <- function(r){
  fold = r['fold']
  index = fold + 2
  r[index] = NA
  return(r)
}


sampleMark <- function(dat, perct){
  # Randomly sample markers from each Chr
  # dat: map file
  # perct: numeric between 0 and 1.
  dat_sub <- dat %>% 
    mutate(markerName = rownames(dat)) %>% 
    group_by(chr) %>% 
    sample_n(perct*length(pos)) %>% 
    as.data.frame()
  rownames(dat_sub) <- dat_sub$markerName
  dat_sub <- dat_sub %>% 
    dplyr::select(-markerName)
}

relatedness_corr <- function(mat1, mat2, plot_draw = TRUE){
  nonzero_mat1 <- which(mat1 != 0, arr.ind = T)
  mat1_coef <- as.vector(mat1[nonzero_mat1])
  mat2_coef <- as.vector(mat2[nonzero_mat1])
  r_corr <- cor(mat1_coef, mat2_coef) # 0.611
  
  cor_G <- data.frame(mat_1 = mat1_coef, mat_2 = mat2_coef)
  g <- ggplot(cor_G, aes(x = mat_1, y = mat_2)) + 
    geom_point() + 
    annotate("text", x=0, y=1.8, label= paste("r = ", r_corr), size = 6)
  if (plot_draw){
    g
  }
  return(r_corr)
}

kfoldValidation_pheno <- function(k_fold, pheno_dat){
  # k_fold: numeric, define the number of fold
  # pheno_dat
  # kinship_mat: matrix, marker-based relateness matrix

  new_ph <- as.data.frame(pheno_dat)
  new_ph$fold <- sample(1:k_fold, size = nrow(new_ph), replace = T)
  new_y <- as.matrix(new_ph[, 1, drop = F]) %*% t(rep(1, k_fold))
  colnames(new_y) = paste0('Y', 1:k_fold)
  new_ph <- cbind(new_ph, new_y)
  
  new_ph_2 <- t(apply(new_ph, 1, missingInFold))
  new_ph_df <- data.frame(new_ph_2)
  new_ph_df$pedigree <- factor(rownames(new_ph_df))
  
  return(new_ph_df)
}

gblup_meta <- function(pheno_dat, kinship_mat, list1, list2){
  # pheno_dat: phenotypic data, vector
  # kinship_mat: matrix, marker-based relationship
  # list1: list, model info
  # list2: list, predictions
  
  varcomps <- sapply(list1, FUN = get_varcomps)
  varcomps <- as.data.frame(t(varcomps))
  varcomps$fold <- as.numeric(gsub('Y', '', rownames(varcomps)))
  
  mu <- colMeans(pheno_dat[, paste('Y',  1:k_fold, sep = '')], na.rm = T)
  mu <- as.data.frame(mu)
  mu$fold <- as.numeric(gsub('Y', '', rownames(mu)))
  
  vc_mu <- merge(varcomps, mu)
  preds <- t(apply(pheno_dat, 1, get_pred, list2))
  
  pred_obs <- merge(pheno_dat[,c('LSM.1', 'fold')],
                    as.data.frame(preds), by = 0)
  
  pred_obs <- merge(pred_obs, vc_mu, by = 'fold')
  pred_obs$pred <- pred_obs$mu + pred_obs$solution
  colnames(pred_obs) <- c('fold', 'Pedigree', 'Observed', 'BLUP', 'BLUP.se', 'z', 'Vg', 'Vresid', 'mu', 'Predicted')
  pred_obs$reliability <- with(pred_obs, 1 - (BLUP.se^2*Vresid)/(Vg * (1+mean(diag(kinship_mat))))) 
  pred_obs$h2 <- with(pred_obs, Vg/(Vg+Vresid))
  return(pred_obs)
}
validationAccuracy <- function(dat){
  # dat: data frame with observed phenotypic values and predicted
  # it contains results from K-fold validation
  
  pred_sum <- dat %>%
    group_by(fold) %>%
    mutate(PredCor = cor(Observed, Predicted),
           RankCor = cor(Observed, Predicted, method = 'spearman'),
           MSE = mean(Observed - Predicted)**2,
           Bias = coef(lm(Observed ~ Predicted))[2]) %>% 
    ungroup() %>% 
    dplyr::select(PredCor, RankCor, MSE, Bias)
  
  return(pred_sum)
}

# Genetic distance function from Sean
gDist=function(gen){
  samples=colnames(gen)
  D=matrix(NA,ncol(gen),ncol(gen))
  gen[gen==5]=NA					
  if (class(gen)=='data.frame'){
    ngen=as.matrix(gen)
    ngen=apply(ngen,2,as.numeric)
    gen=ngen
  }
  for (i in 1:ncol(gen)){
    vi=gen[,i]
    gi=abs(vi-gen)/2
    pi=apply(gi,2,sum,na.rm=TRUE)/(nrow(gi)-apply(is.na(gi),2,sum))	#Want a vect of subs per site avg'd #Remember this is going to be 2 alleles values (2,1,0) alleles diff. need to divide by two.
    D[,i]=pi
    D[i,i]=NA
  }
  row.names(D)=samples
  colnames(D)=samples
  D
}

# From Sean's Script
plotEigenPed=function(gen){
  #cat('pca read.table eigenvec file with id, fam and at least the first 8 eigenvec eg gcta and title for plots. Probably want a pdf open.')
  # if ('ape' %in% installed.packages() == FALSE){
  #   install.packages(c('ape','picante'),quiet=T,repos='http://lib.stat.cmu.edu/R/CRAN')}
  library(ape)
  G=gDist(gen)
  gca=pcoa(as.dist(G))
  biplot(gca,main='PCA Plot of all samples')
  ped=data.frame('Type'=NA,'Pedigree'=colnames(gen),'Phen'=NA)
  for (i in 1:nrow(ped)){
    ped[i,]$Type=substr(as.character(ped$Pedigree[i]),1,3)}	#Writes first three chars of Ped col to get MCT type.
  common=unique(ped$Type)[summary(as.factor(ped$Type)) > 1]
  ped=ped[ped$Type %in% common,]
  gen=gen[,ped$Pedigree]
  G=gDist(gen)
  len=length(unique(ped$Type))
  gca=pcoa(as.dist(G))
  pca=cbind(ped[,1:2],gca$vectors)
  pca=pca[,1:10]
  colnames(pca)=c('fam','id','pca1','pca2','pca3','pca4','pca5','pca6','pca7','pca8')
  #pca=pca[pca[,1] %in% names(table(pca[,1])[table(pca[,1])>2]),]
  len=length(unique(pca[,1]))
  if (len > 30){
    p=palette(colors(distinct=T))[seq(5,5+(len*3),3)]     #Randomly assign colours. try and jump by three so there not to similar.
  }else{ p=palette(c('red','yellow2','orange','blue','green','brown','cyan','magenta','cadetblue','chartreuse4',
                     'purple','violet','seagreen','aquamarine2','steelblue','pink','lightblue','rosybrown','greenyellow','gray',
                     'black','slategray','plum','wheat2','moccasin','hotpink','gold','tan1','violetred','honeydew'))}
  plot(pca$pca1,pca$pca2,col=as.factor(pca[,1]),xlab='pca1',ylab='pca2', lwd=2, main='PCA Plot of main MCTs')
  legend( x='bottomleft', legend=unique(pca[,1]), col=unique(as.factor(pca[,1])), pch=1, bty='n', cex=0.5)
  #plot(pca$pca2,pca$pca3,col=as.factor(pca[,1]),xlab='pca2',ylab='pca3', lwd=2,main=ttl)
  #legend( x='bottomleft', legend=unique(pca[,1]), col=unique(as.factor(pca[,1])),pch=1, bty='n' )
}


gen2AddiFast=function(gen){
  #gena=t(apply(gen[,9:ncol(gen)],1,function(x)sapply(apply(data.frame(x),2,function(y)strsplit(y,''))[[1]],function(z)sum(x[[1]]==strsplit(z,'')[[1]]))))
  #Almost works! 
  #tmp=apply(gen[,9:ncol(gen)],1,function(x)sapply(apply(data.frame(x),2,function(y) if(length(y) > 0)strsplit(y,'')[[1]]),function(z)sum(z==strsplit(x[1],'')[[1]])))
  #A little faster. Would be ideal to do in all apply functions. 
  #tmp=apply(gen[,9:ncol(gen)],1,function(x) unique(unlist(strsplit(unique(x[!is.na(x)]),''))))
  #print(class(gen))
  #print(head(gen[,c(1,2,3,4,9,10,11,12)]))
  allele=apply(gen[,9:ncol(gen)],1,function(x) unique(unlist(strsplit(unique(x[!is.na(x)]),''))))
  ###For some reason some times allele was coming out as matrix rather than alleles?
  if(class(allele)=='matrix') {
    nms=colnames(allele)
    allele=split(allele, rep(1:ncol(allele), each = nrow(allele)))
    names(allele)=nms}
  allele=allele[!sapply(allele,function(x) is.null(x)|'*' %in% x)]
  gena=data.frame(matrix(NA,nrow=length(allele[!sapply(allele,is.null)]),ncol=ncol(gen[,9:ncol(gen)])))
  #print(head(names(allele)))
  rownames(gena)=names(allele)
  colnames(gena)=colnames(gen[,9:ncol(gen)])
  for (i in names(allele)){
    calls=names(sort(table(factor(as.character(gen[i,9:ncol(gen)])))))
    maj=unlist(strsplit(calls[length(calls)],''))[1]
    gena[i,]=apply(gen[i,9:ncol(gen)],2,function(x) sum(strsplit(x,'')[[1]]==maj))
  }
  gena
}
