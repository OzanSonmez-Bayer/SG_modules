# Functions to calculate A matrix

source('/repos/ITG_Data/Code/Shared/deepped.R')


Amat_2gen <- function(dat) {
  
  if (!'P1' %in% colnames(dat)){
    print("parentage pedigree doesn't exit")
    dat <- dat %>% 
      separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
    p1 <- unlist(lapply(dat$P1, FUN = pedSingle))
    p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    dat$P1 <- p1
    dat$P2[dat$P2=='']=NA
    p2 <- unlist(lapply(dat$P2, FUN = pedSingle))
    p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    dat$P2 <- p2
  }
  raw_dat_sub <- dat %>% 
    mutate(P1 = as.character(P1), 
           P2 = as.character(P2),
           trait = OBSRVTN_REF_CD, 
           TRAIT_VALUE = as.numeric(TRAIT_VALUE), 
           GROWSEASON = as.character(GROWSEASON),
           pedigree = as.character(PEDIGREE_NAME))
  raw_dat_sub[is.na(raw_dat_sub$P1), 'P2'] <- NA
  raw_dat_sub[is.na(raw_dat_sub$P2), 'P1'] <- NA
  
  
  ped_file <- raw_dat_sub %>% 
    dplyr::select(PEDIGREE_NAME, P1, P2) %>% 
    distinct_all()
  
  # Two generation A matrix
  
  ped_file$Selfing <- rep(0, rep = nrow(ped_file))
  A_2gen <- asreml.Ainverse(ped_file, fgen = c('Selfing',5))
  A_inv_2gen <- A_2gen$ginv
  A_inv_sparse_2gen <- asreml.sparse2mat(A_inv_2gen)
  Additive_2gen <- round(solve(A_inv_sparse_2gen),2) 
  
  
  colnames(Additive_2gen) <- A_2gen$pedigree$PEDIGREE_NAME
  rownames(Additive_2gen) <- A_2gen$pedigree$PEDIGREE_NAME
  
  return(Additive_2gen)
}




Amat_5gen <- function(dat, n_gen) {
  
  if (!'P1' %in% colnames(dat)){
    print("parentage pedigree doesn't exit")
    dat <- dat %>% 
      separate(ORIGIN, c('P1','P2'), sep = '[\\+]', remove = FALSE, extra = 'merge', fill = 'left')
    p1 <- unlist(lapply(dat$P1, FUN = pedSingle))
    p1 <- sapply(strsplit(p1,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p1 <- sapply(strsplit(p1,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p1 <- sapply(strsplit(p1,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    dat$P1 <- p1
    dat$P2[dat$P2=='']=NA
    p2 <- unlist(lapply(dat$P2, FUN = pedSingle))
    p2 <- sapply(strsplit(p2,"^[A-Z][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p2 <- sapply(strsplit(p2,"^[A-Z][0-9]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    p2 <- sapply(strsplit(p2,"^[0-9][A-Z]_"),function(x)if(length(x)>1)x[[2]] else x[[1]])
    dat$P2 <- p2
  }
  raw_dat_sub <- dat %>% 
    mutate(P1 = as.character(P1), 
           P2 = as.character(P2),
           trait = OBSRVTN_REF_CD, 
           TRAIT_VALUE = as.numeric(TRAIT_VALUE), 
           GROWSEASON = as.character(GROWSEASON),
           pedigree = as.character(PEDIGREE_NAME))
  raw_dat_sub[is.na(raw_dat_sub$P1), 'P2'] <- NA
  raw_dat_sub[is.na(raw_dat_sub$P2), 'P1'] <- NA
  
  
  # ped_file <- raw_dat_sub %>% 
  #   dplyr::select(PEDIGREE_NAME, P1, P2) %>% 
  #   distinct_all()
  
  new_token(fingerprint = F)
  gen5 <- make_ped_file(raw_dat_sub, crop_ids = CropIDs, n_gen = n_gen)
  A_matrix5 <- asreml.Ainverse(gen5$ped, fgen = list("inbreeding", n_gen))
  
  A_inv_5gen <- A_matrix5$ginv
  A_inv_sparse_5gen <- asreml.sparse2mat(A_inv_5gen)
  Additive_5gen <- round(solve(A_inv_sparse_5gen),2) 
  
  # A_matrix5_1 <- left_join(A_matrix5$pedigree, gen5$distinct_lines, by  = 'ID')
  
  # A_matrix5_1 <- inner_join(A_matrix5$pedigree, gen5$distinct_lines, by  = 'ID')
  colnames(Additive_5gen) <- A_matrix5$pedigree$ID
  rownames(Additive_5gen) <- A_matrix5$pedigree$ID
  # print('Additive_5gen')
  # print(dim(Additive_5gen))
  # 
  # Additive_5gen_sub <- Additive_5gen[rownames(Additive_5gen) %in% A_matrix5_1$ID, colnames(Additive_5gen) %in% A_matrix5_1$ID]
  # A_matrix5_1 <- A_matrix5_1[match(rownames(Additive_5gen_sub), A_matrix5_1$ID), ]
  
  ped_name <- CropIDs %>% 
    dplyr::select(M.GERMPLASM.X_ID, M.GERMPLASM.PEDIGREE) %>% 
    dplyr::rename(ID = M.GERMPLASM.X_ID) %>% 
    unique()
  
  ped_name <- inner_join(A_matrix5$pedigree, ped_name, by = 'ID')
  A_matrix5_1 <- Additive_5gen[rownames(Additive_5gen) %in% ped_name$ID, colnames(Additive_5gen) %in% ped_name$ID]
  A_matrix5_2 <- A_matrix5_1[match(rownames(A_matrix5_1), ped_name$ID), ]
  
  
  rownames(A_matrix5_2) <- ped_name$M.GERMPLASM.PEDIGREE
  colnames(A_matrix5_2) <- ped_name$M.GERMPLASM.PEDIGREE
  
  
  return(A_matrix5_2)
}