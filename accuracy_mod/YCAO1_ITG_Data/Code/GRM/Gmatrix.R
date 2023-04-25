# Valildate Genomic Data

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

# source('Code/Shared/BLUPFunctions.R')
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
  dplyr::rename(chr = Chr, 
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

trt <- 'FBRIX'
pheno_df <- pheno %>% 
  filter(trait == trt) %>% 
  mutate(dEBV = predicted.value/reliability)

length(intersect(as.character(pheno_df$PEDIGREE_NAME), as.character(geno$Pedigree))) # 496 inbreds have been genotyped

# Match Geno with pheno

geno_sub <- geno %>% 
  filter(as.character(Pedigree) %in% as.character(pheno_df$PEDIGREE_NAME))
rownames(geno_sub) <- as.character(geno_sub$Pedigree)
geno_sub <- geno_sub[,-1]
geno_cln <- genClean(geno_sub)
rownames(geno_cln) <- rownames(geno_sub)

pheno_sub <- pheno_df %>% 
  filter(as.character(PEDIGREE_NAME) %in% rownames(geno_sub)) %>% 
  dplyr::select(PEDIGREE_NAME, dEBV)
pheno_sub <- pheno_sub[!duplicated(pheno_sub$PEDIGREE_NAME),]
rownames(pheno_sub) <- pheno_sub$PEDIGREE_NAME
pheno_sub <- pheno_sub %>% 
  dplyr::select(-PEDIGREE_NAME)


############## Kinship Matrix ##########################

# Create gpData Object

gp <- create.gpData(geno = geno_cln, pheno = pheno_sub, pedigree = NULL, map = NULL)
summary(gp)

# Recode loci and impute missing genotypes
gp_num <- codeGeno(gp, 
                   label.heter = 'alleleCoding', 
                   maf = 0.05,
                   nmiss = 0.25,
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


###################### Different Ways to Calcualte Unweighted G Matrix #################

# Calculate by Equation
M <- gp_num$geno
p <- apply(M, 2, mean)/2
P <- matrix(rep(p*2, nrow(M)), ncol = ncol(M), nrow = nrow(M), byrow = TRUE)
rownames(P) <- rownames(M)
colnames(P) <- colnames(M)
Z <- M - P
q <- 1 - p
sum2pq <- 2*sum(p*q)
print(sum2pq)
G <- (Z %*% t(Z))/sum2pq

G <- G[rownames(G_base), colnames(G_base)]

cor(diag(G), diag(G_base))
cor(as.vector(G), as.vector(G_base))
cor(G[lower.tri(G, F)],G_base[lower.tri(G_base, F)]) # 0.8 2.8

# Calculate by snpReady

G_snpReady_V <- snpReady::G.matrix(gp_num$geno, method = 'VanRaden', format = 'wide')
G_snpReady_V <- G_snpReady_V$Ga
hist(diag(G_snpReady_V))

cor(diag(G_snpReady_V), diag(G_base))
cor(as.vector(G_snpReady_V), as.vector(G_base))

# Calculate by snpReady with Yang's method

G_snpReady_Y <- snpReady::G.matrix(gp_num$geno, method = 'UAR', format = 'wide')
hist(diag(G_snpReady_Y))

cor(diag(G_snpReady_Y), diag(G_base)) # 0.8642
cor(as.vector(G_snpReady_Y), as.vector(G_base)) # 0.9858

# Calculate by AGHmatrix with Yang
G_AGHmatrix_Y <- AGHmatrix::Gmatrix(gp_num$geno, method = 'Yang')
hist(diag(G_AGHmatrix_Y))

cor(diag(G_snpReady_Y), diag(G_AGHmatrix_Y)) # 1
cor(as.vector(G_snpReady_Y), as.vector(G_AGHmatrix_Y)) # 1


# Calculate by AGHmatrix with VanRaden

G_AGHmatrix_V <- AGHmatrix::Gmatrix(gp_num$geno, method = 'VanRaden')
hist(diag(G_AGHmatrix_V))

cor(diag(G_base), diag(G_AGHmatrix_V)) # 1
cor(as.vector(G_base), as.vector(G_AGHmatrix_V)) # 0.99

# Calculate by rrBLUP with {-1, 0, 1}
G_rrBLUP <- rrBLUP::A.mat((gp_num$geno - 1))
hist(diag(G_rrBLUP))
range(diag(G_rrBLUP)) # 0.80 2.89

cor(diag(G_base), diag(G_rrBLUP)) # 1
cor(as.vector(G_base), as.vector(G_rrBLUP)) # 0.99


#################################### Use Full Taqman Data ##################################

procTom_taqman_2[1:5,1:5]

taqman_full <- procTom_taqman_2 %>% 
  dplyr::select(-Inventory_SampleID) %>% 
  unique()
rownames(taqman_full) <- taqman_full$Pedigree
taqman_full <- taqman_full %>% 
  dplyr::select(-Pedigree)

# Use snpReady to QA/QC
taqman_cln <- genClean(taqman_full)
rownames(taqman_cln) <- rownames(taqman_full)
# taqman_cln <- raw.data(data = as.matrix(taqman_cln), 
#                        frame = 'wide',
#                        base = TRUE, 
#                        sweep.sample = 0.5,
#                        call.rate = 0.95,
#                        maf = 0.05, 
#                        imput = FALSE, 
#                        imput.type = 'wright',
#                        outfile = '012')
# taqman_cln_1 <- taqman_cln$M.clean # 2581 by 253
# 
# # results from snpReady and rrBLUP are different, Why? Which one? 
# G_V <- G.matrix(M = taqman_cln_1, method = 'VanRaden', format = 'wide')$Ga
# hist(diag(G_V))
# range(diag(G_V)) # 0.3788, 1.777

taqman_full_1 <- gen2Additive(taqman_full)
rownames(taqman_full_1) <- rownames(taqman_full)

G_rrBLUP <- rrBLUP::A.mat((taqman_full_1 - 1))
hist(diag(G_rrBLUP))
range(diag(G_rrBLUP)) 


G_Y <- G.matrix(M = taqman_cln_1, method = 'UAR', format = 'wide')
hist(diag(G_Y))
range(diag(G_Y))# 0.77 3.24 center around 2

# Calculate G by equation from VanRaden {0, 1, 2}
M_full <- taqman_cln_1
p_full <- apply(M_full, 2, mean)/2
P_full <- matrix(rep(p_full*2, nrow(M_full)), ncol = ncol(M_full), nrow = nrow(M_full), byrow = TRUE)
rownames(P_full) <- rownames(M_full)
colnames(P_full) <- colnames(M_full)
Z_full <- M_full - P_full
q_full <- 1 - p_full
sum2pq_full <- 2*sum(p_full*q_full)
print(sum2pq_full) # 92.4446, close to the scaling factor calculated by 496 pedigrees, 94.57874
G_full <- (Z_full %*% t(Z_full))/sum2pq_full

hist(diag(G_full))
range(diag(G_full)) # 0.719, 3.373

############################## Check FP Data ##########################
AWS.Authorization('ycao1')
path <- 'ycao1/ITG_data/ProcTomato/'
bucketname <- 'genome-analytics-perm-space'
infile <- paste0('s3://', bucketname, '/', path, 'Geno/', 'Infinium.csv')
csvfile <- get_object(infile)
csvfile_1 <- rawToChar(csvfile)
con <- textConnection(csvfile_1)
FP_geno <- read.csv(con)
close(con)

# FP_geno <- read.csv('fingerprint_data.csv')
FP_map <- FP_geno[,1:3]
FP_geno_1 <- FP_geno[, -c(1:9)]
ped_na <- colnames(FP_geno_1)
FP_geno_1 <- t(FP_geno_1)
colnames(FP_geno_1) <- as.vector(FP_geno[,1])

# FP_geno_cln <- genClean(FP_geno_1)
# FP_geno_cln_1 <- raw.data(data = as.matrix(FP_geno_cln),
#                                         frame = 'wide',
#                                         base = TRUE,
#                                         sweep.sample = 0.5,
#                                         call.rate = 0.95,
#                                         maf = 0.05,
#                                         imput = TRUE,
#                                         imput.type = 'wright',
#                                         outfile = '012')
# 
# FP_sub_cln_ready <- FP_sub_cln_1$M.clean

FP_miss_mrk <- apply(FP_geno_1, MARGIN = 1, FUN = function(x){sum(is.na(x))/ncol(FP_geno_1)})
FP_miss_ind <- apply(FP_geno_1, MARGIN = 2, FUN = function(x){sum(is.na(x))/nrow(FP_geno_1)})

# PSQ_FP <- FP_geno[substr(rownames(FP_geno), start = 1, stop = 3) == 'PSQ', ]

taqman_cln_1 <- taqman_cln
taqman_ped <- gsub('\\-', '.', rownames(taqman_cln_1))
taqman_ped <- gsub('\\_', '.', taqman_ped)

length(intersect(rownames(FP_geno_1), taqman_ped))

shared_ped <- intersect(rownames(FP_geno_1), taqman_ped)
# [1] "PSQ.XL16.1076"  "PSQXL15.1058"   "PSQ.9Z15.9125"  "PSQXL12.1016"   "PSQ.9Z06162"    "PSQ.24.2082"    "PSQ.XL16.1204"  "PSQ.9Z09029"   
# [9] "PSQXL14.1032"   "PSQ9Z11.9049"   "PSQ.176.DAMLER" "PSQ.9Z15.1608V" "PSQ.XL16.1203"  "PSQXL15.1061"   "PSQXL14.1046"   "PSQ.XL17.1087" 
# [17] "PSQ.9Z05086"    "PSQXL14.1037"   "PSQ.9Z15.3843V" "PSQXL14.1041"   "PSQ.9Z18.9175"  "PSQXL14.1045"   "PSQ.XL18.1093"  "PSQXL14.1044"  

FP_geno_1[which(FP_geno_1 == '??')] <- NA
# FP_sub_cln <- gen2Additive(FP_geno_1)

FP_sub <- FP_geno_1[rownames(FP_geno_1) %in% shared_ped,]
FP_sub_cln <- genClean(FP_geno_1)
rownames(FP_sub_cln) <- rownames(FP_geno_1)

# FP_sub_cln_1 <- raw.data(data = as.matrix(FP_sub_cln), 
#                        frame = 'wide',
#                        base = TRUE, 
#                        sweep.sample = 0.5,
#                        call.rate = 0.95,
#                        maf = 0.05, 
#                        imput = TRUE, 
#                        imput.type = 'wright',
#                        outfile = '012')
# FP_sub_cln_ready <- FP_sub_cln_1$M.clean # 4582 5620

# Taqman Subset

# taqman_cln <- taqman_full_1
rownames(taqman_cln) <- gsub('\\_', '.', rownames(taqman_cln))
rownames(taqman_cln) <- gsub('\\-', '.', rownames(taqman_cln))

rownames(taqman_cln_1) <- taqman_ped

Taqman_sub <- taqman_cln[rownames(taqman_cln) %in% shared_ped, ]
FP_sub_cln_ready_1 <- FP_sub_cln[rownames(FP_sub_cln) %in% shared_ped, ]

Taqman_sub <- Taqman_sub[rownames(FP_sub_cln_ready_1),] # 24, 253




############################# FP vs Taqman ############################

# # PCA
# taqman_t <- t(Taqman_sub)
# GD_taqman <- gDist(taqman_t)
# 
# plotEigenPed(taqman_t)
# 
# 
# FP_t <- t(FP_sub_cln_ready[, colnames(FP_sub_cln_ready) %in% colnames(Taqman_sub)])
# GD_FP <- gDist(FP_t)
# plotEigenPed(GD_FP)
# 
# 
# # Compare GD 
# GD_taqman_lower_tri <- GD_taqman[lower.tri(GD_taqman, diag = F)]
# GD_FP_lower_tri <- GD_FP[lower.tri(GD_FP, diag = F)]
# 
# cor(GD_taqman_lower_tri, GD_FP_lower_tri)
# 
# GD_diff <- abs(GD_taqman - GD_FP)                       
# 
# GD_diff_tri <- GD_diff[lower.tri(GD_diff, diag = F)]
# hist(GD_diff_tri)
# range(GD_diff_tri)  # 0.0002630646 0.1650720929
# 
# # Check lines with more than 15% misalignment 
# which(GD_diff >= 0.15, arr.ind = T)

FP_sub_cln_ready_sub <- FP_sub_cln_ready_1[,colnames(FP_sub_cln_ready_1)[colnames(FP_sub_cln_ready_1) %in% colnames(Taqman_sub)]]
Taqman_sub_sub <- Taqman_sub[, colnames(Taqman_sub)[colnames(Taqman_sub) %in% colnames(FP_sub_cln_ready_sub)]]

rownames(Taqman_sub_sub) <- paste(rownames(Taqman_sub_sub), '_Taq', sep = "")

FP_sub_cln_ready_sub <- FP_sub_cln_ready_sub[, colnames(Taqman_sub_sub)]

FP_taq <- rbind(FP_sub_cln_ready_sub, Taqman_sub_sub)

FP_taq_cln <- gen2Additive(FP_taq)
# FP_taq_cln <- genClean(FP_taq)
rownames(FP_taq_cln) <- rownames(FP_taq)


# FP_taq_012 <- raw.data(as.matrix(FP_taq_cln), 
#                        frame = 'wide',
#                        base = TRUE, 
#                        sweep.sample = 0.5,
#                        call.rate = 1,
#                        maf = 0.01, 
#                        imput = FALSE, 
#                        outfile = '012')


# FP_taq_012 <- gen2AddiFast(FP_taq_cln)


# FP_taq_dist <- gDist(t(FP_taq_012$M.clean))
FP_taq_dist <- gDist(t(FP_taq_cln))

range(diag(FP_taq_dist[1:24, 25:48]))
hist(diag(FP_taq_dist[1:24, 25:48]), xlab = 'Genetic Distance on Diagonal', main = '')
range(FP_taq_dist[1:24, 25:48][lower.tri(FP_taq_dist[1:24, 25:48], diag = F)])

which(diag(FP_taq_dist[1:24, 25:48]) >= 0.02)

# Two samples may have genotyping errors, share with Cristian and Tyler, worth double check. 

# "PSQ9Z11.9049_Taq"
# 
# "PSQXL14.1041_Taq"



# export coefficients on diagonal
dig_df <- data.frame(pedigree = rownames(FP_taq_dist)[1:24], distance = diag(FP_taq_dist[1:24, 25:48]))

########################### FP from Deep Pedigrees ###########################
deep_ped <- gsub('\\-', '.', unique(pheno$PEDIGREE_NAME))
deep_ped<- gsub('\\_', '.', deep_ped)
length(intersect(deep_ped, rownames(FP_geno_1)))

deep_ped_FP <- FP_geno_1[which(rownames(FP_geno_1) %in% deep_ped), ]

hist(FP_miss_mrk[which(names(FP_miss_mrk) %in% deep_ped)])

sum(FP_miss_mrk[which(names(FP_miss_mrk) %in% deep_ped)] < 0.20) # 249 lines


########################## Compare GRMs for 249 Ancestor Lines ################

# FP_GRM <- G.matrix(FP_sub_cln_ready, method = 'VanRaden', format = 'wide')$Ga
# Taqman_GRM <- 
FP_geno_1[which(FP_geno_1 == '??', arr.ind = T)] <- NA
FP_geno_full_1 <- genClean(FP_geno_1)
rownames(FP_geno_full_1) <- rownames(FP_geno_1)
FP_geno_full_2 <- raw.data(data = as.matrix(FP_geno_full_1), 
                         frame = 'wide',
                         base = TRUE, 
                         sweep.sample = 0.5,
                         call.rate = 0.95,
                         maf = 0.05, 
                         imput = TRUE, 
                         imput.type = 'wright',
                         outfile = '012')
dim(FP_geno_full_2$M.clean)

FP.ready <- FP_geno_full_2$M.clean

FP.ready_sub <- FP.ready[which(rownames(FP.ready) %in% deep_ped),]

mono_mrk <- names(apply(FP.ready_sub, MARGIN = 2, 
                        FUN = function(x){length(unique(x))})[apply(FP.ready_sub, MARGIN = 2, 
                                                                    FUN = function(x){length(unique(x))}) == 1])
FP.ready_sub <- FP.ready_sub[, which(!colnames(FP.ready_sub) %in% mono_mrk)]

FP_G_snpReady_Y <- snpReady::G.matrix(FP.ready_sub, method = 'UAR', format = 'wide', plot = F)


FP_G_snpReady_V <- snpReady::G.matrix(FP.ready_sub, method = 'VanRaden', format = 'wide', plot = F)
hist(diag(FP_G_snpReady_V$Ga))
range(diag(FP_G_snpReady_V$Ga)) # 0.59 2.26


realizedPD <- nearPD(FP_G_snpReady_V$Ga, keepDiag = T)
G_base <- matrix(realizedPD[[1]]@x, nrow = realizedPD[[1]]@Dim[1])
G_base <- G_base + diag(0.01, nrow = nrow(G_base))
attr(G_base, 'dimnames') <- realizedPD[[1]]@Dimnames
class(G_base) <- 'relationshipMatrix'
str(G_base)
summary(G_base)

plot(G_base, main = 'Unweighted G Matrix')


# Additive Matrix
ped_mat$ped <- gsub('\\-', '.', ped_mat$Pedigree)
ped_mat$ped <- gsub('\\_', '.', ped_mat$ped)
ped_mat_sub <- ped_mat %>% 
  filter(ped %in% rownames(FP_G_snpReady_Y))

dim(ped_mat_sub)

# ped_file <- create.pedigree(ped_mat_sub$ped, ped_mat_sub$par1, ped_mat_sub$par2, ped_mat_sub$gen2go, ped_mat_sub$inbreeding)
# ped_file_df <- as.data.frame(ped_file)
# colnames(ped_file_df)[4] <- 'inbreeding'

ped_mat_sub_1 <- ped_mat_sub %>% 
  dplyr::select(Pedigree, par1, par2, inbreeding, gen2go) %>% 
  rename( 
    PARENT_FEMALE = par1,
    PARENT_MALE = par2)

A_matrix5 <- asreml.Ainverse(ped_mat_sub_1, fgen = list("inbreeding", 5))


# ped_mat_sub_1 <- ped_mat %>% 
#   dplyr::select(Pedigree, par1, par2, inbreeding) %>% 
#   rename( 
#     PARENT_FEMALE = par1,
#     PARENT_MALE = par2)
# 
# A_matrix5 <- asreml.Ainverse(ped_mat_sub_1, fgen = list("inbreeding", 5))



# A_matrix5 <- asreml.Ainverse(ped_file_df, fgen = list("inbreeding", 5))
A_inv <- A_matrix5$ginv
A_inv_sparse <- asreml.sparse2mat(A_inv)
Additive <- round(solve(A_inv_sparse),2)

id_list <- gsub('\\_', '.', A_matrix5$pedigree$Pedigree)
id_list <- gsub('\\-', '.', id_list)
# colnames(Additive) <- A_matrix5$pedigree$ID
# rownames(Additive) <- A_matrix5$pedigree$ID

colnames(Additive) <- id_list
rownames(Additive) <- id_list

Additive_sub <- Additive[rownames(Additive) %in% rownames(FP_G_snpReady_Y), colnames(Additive) %in% rownames(FP_G_snpReady_Y)]
dim(Additive_sub)

hist(diag(Additive_sub))

attr(Additive_sub, 'dimnames') <- list(rownames(Additive_sub))
class(Additive_sub) <- 'relationshipMatrix'
plot(Additive_sub)

boxplot(as.vector(FP_G_snpReady_V$Ga) ~ as.vector(Additive_sub), xlab = 'A', ylab = 'G')

Atab <- write.relationshipMatrix(Additive_sub, sorting = 'ASReml', type = 'none')
colnames(Atab)[3] <- 'A'
head(Atab)


Gtab <- write.relationshipMatrix(FP_G_snpReady_V$Ga, sorting = 'ASReml', type = 'none')
colnames(Gtab)[3] <- 'G'
head(Gtab)

AG <- merge(Atab, Gtab, all = T)
AG[is.na(AG)] <- 0

AG$diff_coef <- AG$A - AG$G
AG <- AG[order(AG$diff), ]
head(AG)

# A_matrix5$pedigree[A_matrix5$pedigree$ID %in% rownames(Additive_sub)[c(154, 69)],]
# snp_check <- FP.ready_sub[which(rownames(FP.ready_sub) %in% A_matrix5$pedigree[A_matrix5$pedigree$ID %in% rownames(Additive_sub)[c(154, 69)],'ID']), ]
# sum(snp_check[1, ] == snp_check[2, ], na.rm = T)/ncol(snp_check)

A_matrix5$pedigree[id_list %in% rownames(Additive_sub)[c(154, 69)],]
snp_check <- FP.ready_sub[which(rownames(FP.ready_sub) %in% id_list[id_list %in% rownames(Additive_sub)[c(154, 69)]]), ]
sum(snp_check[1, ] == snp_check[2, ], na.rm = T)/ncol(snp_check)


tail(AG)

A_matrix5$pedigree[id_list %in% rownames(Additive_sub)[c(131, 130)],]
snp_check <- FP.ready_sub[which(rownames(FP.ready_sub) %in% id_list[id_list %in% rownames(Additive_sub)[c(131, 130)]]), ]
sum(snp_check[1, ] == snp_check[2, ], na.rm = T)/ncol(snp_check)

# A_matrix5$pedigree[A_matrix5$pedigree$ID %in% rownames(Additive_sub)[c(89, 86)],]
# snp_check <- FP.ready_sub[which(rownames(FP.ready_sub) %in% A_matrix5$pedigree[A_matrix5$pedigree$ID %in% rownames(Additive_sub)[c(89, 86)],'ID']), ]
# sum(snp_check[1, ] == snp_check[2, ], na.rm = T)/ncol(snp_check)

nonzeor_A <- which(Additive_sub != 0, arr.ind = T)
G_base_coef <- as.vector(FP_G_snpReady_V$Ga[nonzeor_A])
A_mat_coef <- as.vector(Additive_sub[nonzeor_A])
cor(G_base_coef, A_mat_coef) # 0.611

cor_G_A <- data.frame(G_coef = G_base_coef, A_coef = A_mat_coef)
ggplot(cor_G_A, aes(x = G_coef, y = A_coef)) + 
  geom_point() + 
  annotate("text", x=0, y=1.8, label= "r = 0.83", size = 6)

############################## Shared Markers Between Taqman and FP #################

# With clean taqman and FP data

shared_mrk <- intersect(colnames(Taqman_sub), colnames(FP_sub_cln_ready))
length(shared_mrk)

# Match column order
FP_compare <- FP_sub_cln_ready_1[,colnames(FP_sub_cln_ready_1) %in% shared_mrk]
Taqman_compare <- Taqman_sub[, colnames(Taqman_sub) %in% shared_mrk]

identical(colnames(FP_compare), colnames(Taqman_compare))


FP_compare <- FP_compare[, shared_mrk]
Taqman_compare <- Taqman_compare[, shared_mrk]

identical(colnames(FP_compare), colnames(Taqman_compare))
identical(rownames(FP_compare), rownames(Taqman_compare))

compare_mat <- FP_compare - Taqman_compare

FP_compare[which(compare_mat != 0, arr.ind = T)]
Taqman_compare[which(compare_mat != 0, arr.ind = T)]

length(FP_compare[which(compare_mat != 0, arr.ind = T)]) # 122 different genotypes from Taqman and FP 

length(which(abs(compare_mat) == 2)) # 63

length(which(abs(compare_mat) == 1)) # 59 heterozygous and homozygous switch


# Correct different code 

FP_sub_cln_ready_sub_1 <- FP_sub_cln_ready_sub
FP_sub_cln_ready_sub_1[which(abs(compare_mat) == 2, arr.ind = T)] <- 2
FP_sub_cln_ready_sub_1[which(abs(compare_mat) == 1, arr.ind = T)] <- 1

Taqman_sub_sub_1 <- Taqman_sub_sub
Taqman_sub_sub_1[which(abs(compare_mat) == 1, arr.ind = T)] <- 1
Taqman_sub_sub_1[which(abs(compare_mat) == 2, arr.ind = T)] <- 2
FP_taq <- rbind(FP_sub_cln_ready_sub_1, Taqman_sub_sub_1)
FP_taq_dist <- gDist(t(FP_taq))
range(diag(FP_taq_dist[1:24, 25:48]))
range(FP_taq_dist[1:24, 25:48][lower.tri(FP_taq_dist[1:24, 25:48], diag = F)])







# Correlation between off-diagonal coefficients

load('Additive_2gen.RData')

# A_2gen vs G

A_2gen_sub <- Additive_2gen[which(rownames(Additive_2gen) %in% rownames(G_rrBLUP)),
                            which(colnames(Additive_2gen) %in% colnames(G_rrBLUP))]
dim(A_2gen_sub)

G_sub <- G_rrBLUP[which(rownames(G_rrBLUP) %in% rownames(A_2gen_sub)),
                  which(colnames(G_rrBLUP) %in% colnames(A_2gen_sub))]

A_2gen_sub <- A_2gen_sub[rownames(G_sub), colnames(G_sub)]


A_lowtri_coef <- A_2gen_sub[lower.tri(A_2gen_sub, diag = F)]
G_lowtri_coef <- G_sub[lower.tri(G_sub, diag = F)]

plot(A_lowtri_coef, G_lowtri_coef)


# A_5gen vs G

load('Additive_5gen.RData')



G_rrBLUP_1 <- G_rrBLUP
# G_rrBLUP_1_sub <- G_rrBLUP_1[rownames(G_rrBLUP_1) %in% unique(pheno_df$PEDIGREE_NAME), 
#                              colnames(G_rrBLUP_1) %in% unique(pheno_df$PEDIGREE_NAME)]
# G_rrBLUP_1 <- G_rrBLUP_1_sub
rownames(G_rrBLUP_1) <- gsub('\\_', '.', rownames(G_rrBLUP_1))
rownames(G_rrBLUP_1) <- gsub('\\-', '.', rownames(G_rrBLUP_1))
colnames(G_rrBLUP_1) <- rownames(G_rrBLUP_1)

A_5gen_sub <- Additive_5gen[which(rownames(Additive_5gen) %in% rownames(G_rrBLUP_1)),
                            which(colnames(Additive_5gen) %in% colnames(G_rrBLUP_1))]
dim(A_5gen_sub)

G_sub <- G_rrBLUP_1[which(rownames(G_rrBLUP_1) %in% rownames(A_5gen_sub)),
                  which(colnames(G_rrBLUP_1) %in% colnames(A_5gen_sub))]

A_5gen_sub <- A_5gen_sub[rownames(G_sub), colnames(G_sub)]


A_lowtri_coef <- A_5gen_sub[lower.tri(A_5gen_sub, diag = F)]
G_lowtri_coef <- G_sub[lower.tri(G_sub, diag = F)]

plot(A_lowtri_coef, G_lowtri_coef)

cor(A_lowtri_coef, G_lowtri_coef)

coef_df <- data.frame(A_coef = A_lowtri_coef, G_coef = G_lowtri_coef)


ggplot(coef_df, aes(x=G_coef, y=A_coef)) + 
  geom_point()


plot(diag(A_5gen_sub), diag(G_sub))


############################# Calculate the Differences in Coef between G and A #########################
pivot_sym_mat <- function(in_mat) {
  tmp <- in_mat
  
  tmp[lower.tri(tmp, diag = F)] <- NA
  
  G_df <- as.data.frame(tmp)
  
  G_df <- G_df %>% 
    mutate(ped1 = rownames(.))
  
  G_pivot <- G_df %>% 
    pivot_longer(-ped1, names_to = "ped2", values_to = "coef") %>% 
    drop_na()
  G_pivot
}


G_pivot <- pivot_sym_mat(G_rrBLUP)

G_pivot$ped1 <-  gsub('\\_', '.', G_pivot$ped1)
G_pivot$ped1 <-  gsub('\\-', '.', G_pivot$ped1)

G_pivot$ped2 <-  gsub('\\_', '.', G_pivot$ped2)
G_pivot$ped2 <-  gsub('\\-', '.', G_pivot$ped2)


load('Additive_5gen.RData')
A_pivot <-  pivot_sym_mat(Additive_5gen)




G_A_compare <- G_pivot %>% 
  inner_join(A_pivot, by = c('ped1', 'ped2'), suffix = c('_G', '_A'))

G_A_compare <- G_A_compare %>% 
  mutate(
         coef_G = round(coef_G, 1), 
         coef_A = round(coef_A, 1),
         coef_diff = abs(coef_G - coef_A),)
head(G_A_compare)


View(G_A_compare[which(G_A_compare$coef_G < 0.5 & G_A_compare$coef_A == 1 ), ])



