# GBLUP for Processing Tomato

library(readxl)
library(tidyverse)
library(snpReady)
library(synbreed)
library(ggplot2)
library(asreml)
library(asremlPlus)
library(aws.s3)
library(jsonlite)
library(Matrix)

source('Code/Shared/BLUPFunctions.R')
source('Code/Shared/genoFunctions.R')
source('Code/Shared/misFunctions.R')
source('Code/Shared/CGS.R')

################# Import Raw Taqman Data #################
procTom_taqman <- read.table('ProcTomato/9Z_TaqManFromCGS_FilteredMarkers.tab', 
                             sep = '\t',
                             fill = TRUE,
                             header = TRUE)
procTom_taqman <- procTom_taqman %>% 
  filter(!is.na(Pedigree))

dim(procTom_taqman)
miss_mrk <- apply(procTom_taqman[,-c(1:2)], MARGIN = 1, FUN = function(x){sum(is.na(x))})
miss_ind <- apply(procTom_taqman[,-c(1:2)], MARGIN = 2, FUN = function(x){sum(is.na(x))})

# remove individuals with 285 missing markers
procTom_taqman_1 <- procTom_taqman[-which(miss_mrk == 287),] # 2607 individuals left
miss_ind <- apply(procTom_taqman_1[,-c(1:2)], MARGIN = 2, FUN = function(x){sum(is.na(x))})
procTom_taqman_2 <- procTom_taqman_1[, -which(miss_ind == 2599)]
geno <- procTom_taqman_2[,-1]

###################### Pull Tomato Map ########################

source('Credential/Credential_sdi.R')
Sys.setenv(CLIENT_SECRET = FINGERPRINT_CLIENT_SECRET)
Sys.setenv(CLIENT_ID = FINGERPRINT_CLIENT_ID)
Sys.setenv(PING_URL = FINGERPRINT_PING_URL)
createPingToken()

availableTaxaJson <- pingSecuredGet(URLencode('https://genetics.ag/genetic-maps/v1/'))

availableTaxa <- fromJSON(content(availableTaxaJson, as = "text", encoding = "UTF-8"))

MapID <- availableTaxa[availableTaxa$description == 'Tomato Consensus Map, June 2014',1]

map <- RetrieveMap(MapID)
head(map)
colnames(map) <- c('Chr', 'MarkerName', 'Pos')


map <- map %>% 
  filter(MarkerName %in% colnames(procTom_taqman_2)[-c(1:2)])
rownames(map) <- map$MarkerName

map_sub <- map %>% 
  dplyr::select(-MarkerName) %>% 
  rename(chr = Chr, 
         pos = Pos)
############### Import GCA from Deep Pedigree #############
 
AWS.Authorization('ycao1')
path <- 'ycao1/ITG_data/ProcTomato/'
bucketname <- 'genome-analytics-perm-space'
infile <- paste0('s3://', bucketname, '/', path, 'BLUP/', 'PGCA_PCM0_2014_2019_deep.csv')
csvfile <- get_object(infile)
csvfile_1 <- rawToChar(csvfile)
con <- textConnection(csvfile_1)
pheno <- read.csv(con)
close(con)

############## Match Geno with Pheno #########################

# FBRIX

trt <- 'MATR'
pheno_df <- pheno %>% 
  filter(trait == trt) %>% 
  mutate(dEBV = predicted.value/reliability)

length(intersect(as.character(pheno_df$PEDIGREE_NAME), as.character(geno$Pedigree))) # 496 inbreds have been genotyped

# Match Geno with pheno

# Use all the available lines to calculate allele frequencies
geno_sub <- geno %>%
  filter(as.character(Pedigree) %in% as.character(pheno_df$PEDIGREE_NAME))
# geno_sub <- geno[-c(21,81),] # remove duplicated lines. 100% missing markers
rownames(geno_sub) <- as.character(geno_sub$Pedigree)
geno_sub <- geno_sub[,-1]
geno_cln <- genClean(geno_sub)
rownames(geno_cln) <- rownames(geno_sub)
# 
# geno_cln <- genClean(gen = geno)
# geno_cln <- geno_cln[-c(21,81),]
# rownames(geno_cln) <- geno$Pedigree[-c(21,81)]



pheno_sub <- pheno_df %>% 
  filter(as.character(PEDIGREE_NAME) %in% rownames(geno_sub)) %>%
  dplyr::select(PEDIGREE_NAME, dEBV)
pheno_sub <- pheno_sub[!duplicated(pheno_sub$PEDIGREE_NAME),]
rownames(pheno_sub) <- pheno_sub$PEDIGREE_NAME
pheno_sub <- pheno_sub %>% 
  dplyr::select(-PEDIGREE_NAME)


############## Kinship Matrix ##########################

# Create gpData Object

gp <- create.gpData(geno = geno_cln, pheno = pheno_sub, pedigree = NULL, map = map_sub)
summary(gp)

# Recode loci and impute missing genotypes
gp_num <- codeGeno(gp, 
                   label.heter = 'alleleCoding', 
                   maf = 0.01,
                   nmiss = 0.05,
                   impute = TRUE, 
                   impute.type = 'random', 
                   verbose = TRUE)
# step 1  : 5 marker(s) removed with > 25 % missing values 
# step 2  : Recoding alleles 
# step 4  : 0 marker(s) removed with maf < 0 
# step 7  : Imputing of missing values 
# step 7d : Random imputing of missing values 
# approximate run time for random imputation  0.099  seconds 
# step 8  : No recoding of alleles necessary after imputation 
# step 9  : 0 marker(s) removed with maf < 0 
# step 10 : No duplicated markers removed 
# End     : 280 marker(s) remain after the check
# 
# Summary of imputation 
# total number of missing values                : 3744 
# number of random imputations                  : 3744 


# Kinship Matrix, G_BASE RanVaden

G <- kin(gp_num, ret = 'realized')
summary(G)
G <- G[rownames(gp_num$pheno), rownames(gp_num$pheno)]


# G <- G_AGHmatrix_Y

realizedPD <- nearPD(G, keepDiag = T)
G_base <- matrix(realizedPD[[1]]@x, nrow = realizedPD[[1]]@Dim[1])
G_base <- G_base + diag(0.01, nrow = nrow(G_base))
attr(G_base, 'dimnames') <- realizedPD[[1]]@Dimnames
class(G_base) <- 'relationshipMatrix'
str(G_base)
summary(G_base)

plot(G_base, main = 'Unweighted G Matrix')

# G inverse
G_inv <- write.relationshipMatrix(G_base, sorting = 'ASReml', type = 'ginv')
head(attr(G_inv, 'rowNames'))
names(G_inv) <- c('row', 'column', 'coefficient')
head(G_inv)


######################### GBLUP with Cross Validation #################

k_fold <- 5
set.seed(20200226)
new_ph <- as.data.frame(gp_num$pheno)
new_ph$fold <- sample(1:k_fold, size = nrow(new_ph), replace = T)
new_y <- as.matrix(new_ph[ , 1, drop = T]) %*%t(rep(1, k_fold))
colnames(new_y) <- paste0('Y', 1:k_fold)
new_ph <- cbind(new_ph, new_y)


new_ph_2 <- t(apply(new_ph, 1, missingInFold))
new_ph_df <- data.frame(new_ph_2)
new_ph_df$pedigree <- factor(rownames(new_ph_df))

# result.list <- list()
# pred.list <- list()
# 
# for (trait in names(new_ph_df)[3: (3 + k_fold-1)]){
#   asr <- asreml(fixed = get(trait) ~ 1, 
#                 random = ~giv(pedigree, var = T), 
#                 ginverse = list(pedigree = G_inv), 
#                 rcov = ~units, 
#                 data = new_ph_df, 
#                 na.method.X = 'include')
#   result.list[[trait]] <- asr
#   pred.list[[trait]] <- summary(asr, all = T)$coef.rand
#   rm(list = 'asr')
# }
# 
# varcomps <- sapply(result.list, FUN = get_varcomps)
# varcomps <- as.data.frame(t(varcomps))
# varcomps$fold <- as.numeric(gsub('Y', '', rownames(varcomps)))
# 
# sigma2_list <- sapply(result.list, FUN = function(x){x$sigma2})
# sigma2_df <- as.data.frame(sigma2_list)
# sigma2_df$fold <- as.numeric(gsub('Y', '', rownames(sigma2_df)))
# 
# mu <- colMeans(new_ph_df[, paste('Y', 1:k_fold, sep = "")], na.rm = T)
# mu <- as.data.frame(mu)
# mu$fold <- as.numeric(gsub('Y', '', rownames(mu)))
# 
# vc_mu <- merge(varcomps, mu)
# vc_mu <- merge(vc_mu, sigma2_df)
# 
# 
# preds <- t(apply(new_ph_df, 1, get_pred, pred.list))
# pred_obs <- merge(new_ph_df[,c('dEBV.1', 'fold')],
#                   as.data.frame(preds), by = 0)
# 
# 
# pred_obs <- merge(pred_obs, vc_mu, by = 'fold')
# pred_obs$pred <- pred_obs$mu + pred_obs$solution
# colnames(pred_obs) <- c('fold', 'Pedigree', 'Observed', 'BLUP', 'BLUP.se', 'z', 'Vg', 'Vresid', 'mu', 'sigma2_val', 'Predicted')
# pred_obs$reliability <- with(pred_obs, 1 - (BLUP.se*sigma2_val)/(Vg * (1+diag(G_base)))) # relibiality is negative?
# # pred_obs$reliability <- with(pred_obs, 1 - BLUP.se^2/Vg) # Mrode's formula?
# 
# pred_obs$h2 <- with(pred_obs, Vg/(Vg + Vresid)) # heritability is very small, double check?
# 
# pred_sum <- pred_obs %>% 
#   group_by(fold) %>% 
#   mutate(PredCor = cor(Observed, Predicted), 
#          RankCor = cor(Observed, Predicted, method = 'spearman')) %>%
#   ungroup() %>% 
#   dplyr::select(PredCor, RankCor)
# colMeans(pred_sum)
# 
# 
# # Test GBLUP with all the pedigrees
# GBLUP <- gpMod(gp_num, model = 'BLUP', kin = G_base, trait = 1)
# summary(GBLUP)
# cor(GBLUP$y, GBLUP$g, method = 'spearman')
# cor(GBLUP$y, (GBLUP$g + 5.782))

# Test with Asreml
GBLUP_asreml <- asreml(fixed = dEBV.1 ~ 1,
                       random = ~giv(pedigree, var = T), 
                       ginverse = list(pedigree = G_inv), 
                       rcov = ~units, 
                       data = new_ph_df, 
                       na.method.X = 'include',
                       workspace = 320e+6, pworkspace = 320e+6)
pred_asreml <- as.data.frame(predict(GBLUP_asreml, classify = 'pedigree')$predictions$pval)
cor(new_ph_df$dEBV.1, pred_asreml$predicted.value)

varcomp_tab <- summary(GBLUP_asreml)$varcomp
h2 <- varcomp_tab$component[1]/sum(varcomp_tab$component)
h2
rnd_ebv <- data.frame(coef(GBLUP_asreml)$random, GBLUP_asreml$vcoeff$random)
names(rnd_ebv) <- c('BLUP', 'BLUP.se')
rnd_ebv <- rnd_ebv[grep('giv', rownames(rnd_ebv)),]
rnd_ebv$pedigree <- gsub('giv(pedigree, var = T)_', '', rownames(rnd_ebv), fixed = TRUE)
PEV <- rnd_ebv[, 2]*GBLUP_asreml$sigma2
reliability <- 1 - PEV/varcomp_tab$component[1] * (1 + diag(G_base))



############# H Matrix with AGHmatrix #########################

library(AGHmatrix)

# Get A matrix 
# Use full pedigree set to get A matrix, then subset it by trait

path <- ('ycao1/ITG_data/ProcTomato/')
bucketname <- 'genome-analytics-perm-space'

AWS.Authorization('ycao1')
# source('Credential/Credential_yc.R')

Crop <- 'Tomato'
# Get a list of files in the path
dfBucket <- get_bucket_df(bucketname, paste0(path, 'Data'))
file_list <- dfBucket$Key[grep('PCM', dfBucket$Key)]
file_PCM0 <- file_list[grep('PCM0', file_list)]

# Combined PCM 0 data
PCM0_dat <- data.frame()
for (i in 1:length(file_PCM0)){
  print(file_PCM0[i])
  infile <- paste0('s3://','genome-analytics-perm-space/', file_PCM0[i])
  csvfile <- get_object(infile)
  csvfile_1 <- rawToChar(csvfile)
  con <- textConnection(csvfile_1)
  temp <- read.csv(con)
  close(con)
  if(grep(file_PCM0[i], pattern = 'PCM0')>=1){
    if (i == 1){
      PCM0_dat <- temp
    }else{
      PCM0_dat <- plyr::rbind.fill(PCM0_dat, temp)
    }
  }
}


raw_dat_sub <- PCM0_dat %>% 
  # filter(OBSRVTN_REF_CD == 'FBRIX') %>%
  mutate(P1 = as.character(P1), 
         P2 = as.character(P2),
         trait = OBSRVTN_REF_CD, 
         TRAIT_VALUE = as.numeric(TRAIT_VALUE), 
         GROWSEASON = as.character(GROWSEASON),
         pedigree = as.character(PEDIGREE_NAME))
raw_dat_sub[is.na(raw_dat_sub$P1), 'P2'] <- NA
raw_dat_sub[is.na(raw_dat_sub$P2), 'P1'] <- NA

# raw_dat_sub$P1[is.na(raw_dat_sub$P1)] <- 0
# raw_dat_sub$P2[is.na(raw_dat_sub$P2)] <- 0
# 
# ped_file <- raw_dat_sub %>% 
#   dplyr::select(pedigree, P1, P2) %>% 
#   mutate(pedigree = as.factor(pedigree),
#          P1 = as.factor(P1),
#          P2 = as.factor(P2)) %>% 
#   distinct_all() %>% 
#   rename(Ind = pedigree, 
#          Par1 = P1,
#          Par2 = P2)
# 

ped_file <- raw_dat_sub %>% 
  dplyr::select(PEDIGREE_NAME, P1, P2) %>% 
  distinct_all()

# Two generation A matrix

ped_file$Selfing <- rep(0, rep = nrow(ped_file))
A_2gen <- asreml.Ainverse(ped_file, fgen = c('Selfing',5))
A_inv_2gen <- A_2gen$ginv
A_inv_sparse_2gen <- asreml.sparse2mat(A_inv_2gen)
Additive_2gen <- round(solve(A_inv_sparse_2gen),2) # 2041 individuals


colnames(Additive_2gen) <- A_2gen$pedigree$PEDIGREE_NAME
rownames(Additive_2gen) <- A_2gen$pedigree$PEDIGREE_NAME

# G matrix 

geno_cln_1 <- raw.data(data = as.matrix(geno_cln),
                                        frame = 'wide',
                                        base = TRUE,
                                        sweep.sample = 0.5,
                                        call.rate = 0.95,
                                        maf = 0.05,
                                        imput = TRUE,
                                        imput.type = 'wright',
                                        outfile = '012')
geno_ready <- geno_cln_1$M.clean

G_AGHmatrix_V <- AGHmatrix::Gmatrix(geno_ready, method = 'VanRaden')
dim(G_AGHmatrix_V)


# H matrix with AGHmatrix
Hmat <- Hmatrix(A = Additive_2gen, G = G_AGHmatrix_V, markers = geno_ready)

# Inverse matrix in Asreml readable format
Hmat_inv <- solve(Hmat)
Hmat_inv_asreml <- formatmatrix(Hmat_inv, return = TRUE, exclude.0 = TRUE)


###################################### Calculate BLUP with H matrix #####################################

# Subset Hmatrix based on trait

# Use FBRIX as an example

Hmat_sub <- Hmat[rownames(Hmat)[rownames(Hmat) %in% rownames(pheno_sub)],
                 colnames(Hmat)[colnames(Hmat) %in% rownames(pheno_sub)]]

# Hmat_sub_inv <- formatmatrix(Hmat_sub, return = TRUE, exclude.0 = TRUE)
Hmat_sub_inv <- write.relationshipMatrix(Hmat_sub, sorting = 'ASReml', type = 'ginv')

head(attr(Hmat_sub_inv, 'rowNames'))
names(Hmat_sub_inv) <- c('row', 'column', 'coefficient')

pheno_sub_1 <- data.frame(PEDIGREE_NAME = rownames(pheno_sub), 
                          pheno_sub)
GBLUP_asreml_H <- asreml(fixed = dEBV ~ 1,
                       random = ~ped(PEDIGREE_NAME, var = T), 
                       ginverse = list(PEDIGREE_NAME = Hmat_sub_inv), 
                       rcov = ~units, 
                       data = pheno_sub_1, 
                       na.method.X = 'include')
pred_asreml <- as.data.frame(predict(GBLUP_asreml_H, classify = 'PEDIGREE_NAME')$predictions$pval)
cor(pheno_sub_1$dEBV, pred_asreml$predicted.value) # -0.02

varcomp_tab <- summary(GBLUP_asreml_H)$varcomp
h2 <- varcomp_tab$component[1]/sum(varcomp_tab$component)
rnd_ebv <- data.frame(coef(GBLUP_asreml)$random, GBLUP_asreml$vcoeff$random)
names(rnd_ebv) <- c('BLUP', 'BLUP.se')
rnd_ebv <- rnd_ebv[grep('giv', rownames(rnd_ebv)),]
rnd_ebv$pedigree <- gsub('giv(pedigree, var = T)_', '', rownames(rnd_ebv), fixed = TRUE)
PEV <- rnd_ebv[, 2]*GBLUP_asreml$sigma2
reliability <- 1 - PEV/varcomp_tab$component[1] * (1 + diag(G_base))

######################## GBLUP with FP data ################

AWS.Authorization('ycao1')
infile <- paste0('s3://', bucketname, '/', path, 'Geno/', 'Infinium.csv')
csvfile <- get_object(infile)
csvfile_1 <- rawToChar(csvfile)
con <- textConnection(csvfile_1)
FP_geno <- read.csv(con)
close(con)

FP_geno_1 <- FP_geno[, -c(1:9)]
ped_na <- colnames(FP_geno_1)
FP_geno_1 <- t(FP_geno_1)
colnames(FP_geno_1) <- as.vector(FP_geno[,1])

FP_geno_1[which(FP_geno_1 == '??', arr.ind = T)] <- NA
FP_cln <- genClean(FP_geno_1)
rownames(FP_cln) <- rownames(FP_geno_1)


FP_cln_1 <- raw.data(data = as.matrix(FP_cln), 
                         frame = 'wide',
                         base = TRUE, 
                         sweep.sample = 0.5,
                         call.rate = 0.95,
                         maf = 0.05, 
                         imput = TRUE, 
                         imput.type = 'wright',
                         outfile = '012')
FP_cln_ready <- FP_cln_1$M.clean 

FP_G <- snpReady::G.matrix(FP_cln_ready, method = 'VanRaden', format = 'wide')
FP_G <- FP_G$Ga

# Subset by pedigrees

pheno_sub$PEDIGREE_NAME <- gsub('\\_', '.', rownames(pheno_sub))
pheno_sub$PEDIGREE_NAME <- gsub('\\-', '.', rownames(pheno_sub))

length(intersect(pheno_sub$PEDIGREE_NAME, rownames(FP_G))) # only 12 genotypes

pheno_df$ped <- gsub('\\_', '.',pheno_df$PEDIGREE_NAME)
pheno_df$ped <- gsub('\\-', '.', pheno_df$PEDIGREE_NAME)

length(intersect(unique(pheno_df$ped), rownames(FP_G))) # 125 pedigrees that have been genotyped in FP


