# Weighted G Matrix
library(snpReady)
library(tidyverse)
library(BGLR)
library(MASS)
library(synbreed)
# library(synbreed)

source('Code/Shared/genoFunctions.R')

# use pepper as an example for the calculation of Kinship Matrix

geno_dat <- read.delim(file = 'Pepper/Pepper_geno.txt')
map_dat <- read.delim(file = 'Pepper/Pepper_map.txt')

load('Pepper/dEBV.Rdata')
pheno_dat <- pred_df %>% 
  dplyr::select(Pedigree, dEBV) %>% 
  mutate(Pedigree = gsub('-', '.', Pedigree))
################################### Clean geno and map data ####################################

pos_chr <- map_dat %>% 
  dplyr::select(X, CHROMOSOME, POSITION) %>% 
  mutate(X = as.character(X),
         CHROMOSOME = as.character(CHROMOSOME),
         POSITION = as.character(POSITION))


geno <- t(geno_dat)
colnames(geno) <- geno[1,]
X <- rownames(geno)
geno <- geno[-1,]

geno_filter <- filterfun(geno)
misind <- geno_filter[[1]][geno_filter[[1]] < 0.4] # remove individuals with more than 40% missing markers
mismrk <- geno_filter[[2]][geno_filter[[2]] < 0.25] # remove markers with more than 25% missing individuals. 

geno_c <- geno[, colnames(geno) %in% names(mismrk)]
geno_c[which(geno_c == "??", arr.ind = TRUE)] <- NA


# convert AA to 0, 1, 2
geno_convt <- gen2Additive(geno_c)
rownames(geno_convt) <- rownames(geno_c)
geno_filter <- filterfun(geno_convt)
misind <- geno_filter[[1]][geno_filter[[1]] < 0.6 ]
geno_cln <- geno_convt[rownames(geno_convt) %in% names(misind),] # 4589 individuals and 3715 markers


# Use snpReady to QC data
# start with geno raw data

geno[which(geno == "??", arr.ind = TRUE)] <- NA
geno_convt <- gen2Additive(geno) # convert nitrogenous to numeric
rownames(geno_convt) <- rownames(geno)

misind <- apply(geno_convt, 1, FUN = function(x) {sum(is.na(x))/length(x)})
mismark <- apply(geno_convt, 2, FUN = function(x) {sum(is.na(x))/length(x)})
hist(misind)
hist(mismark)
# Remove individuals with more than 60% missing values
# Remove markers with more than  75% missing obs
geno_ready <- raw.data(data = as.matrix(geno_convt), 
                       frame = 'wide', base = FALSE, 
                       sweep.sample = 0.7, 
                       call.rate = 0.25, 
                       maf = 0.1, 
                       imput = TRUE, 
                       imput.type = 'wright', 
                       outfile = '012')


# Wright method has 62% percent accuracy. Imputation basedon mean has intermediate accuracy as well. 

Mrk_cln <- geno_ready$M.clean


##################### G Matrix Calculation Onward ##############################


########### Step 1: Allele Frequencies

allele_p <- apply(Mrk_cln, 2, mean)/2
allele_p <- round(allele_p, 3)
allele_df <- data.frame(marker = names(allele_p), p = unname(allele_p))
ggplot(allele_df, aes(x=p)) +
  geom_histogram(bins=20, color="white") +
  ggtitle("Histogram of allele frequencies (p)") +
  xlab("Allele Frequencies") +
  theme_minimal()


########### Step 2: Create G matrix

# Center Marker matrix using the means (2 times the allele frequencies)
allele_mat <- matrix(rep(allele_p*2, nrow(Mrk_cln)), ncol = ncol(Mrk_cln), nrow = nrow(Mrk_cln), byrow = TRUE)
rownames(allele_mat) <- rownames(Mrk_cln)
colnames(allele_mat) <- colnames(Mrk_cln)

# Create Z matrix, centered Marker matrix

Z <- Mrk_cln - allele_mat

########## Step 3: SNP weights

# Generic G matrix  calculation; G = Z*D*Z_transpose. D is the scaling factor, SNP weights

### Method 1: W = 2*p*(1-p)
# markers contribute to genomic relatedness proportional to the reciprocal of their expected variance

D_base <- diag(1 / ncol(Mrk_cln) * (2 * allele_p * (1 - allele_p)))
colnames(D_base) <- colnames(Mrk_cln)
rownames(D_base) <- colnames(Mrk_cln)

G_base <- Z %*% D_base %*% t(Z)

G_base_sub <- G_base[which(rownames(G_base) %in% as.character(pheno_dat$Pedigree)),
                     which(colnames(G_base) %in% as.character(pheno_dat$Pedigree))]

### Method 2: W = 2*p*(1-p)*mu^2, mu is the allele effect from Bayesian RR and BL
# This method could be trait specific by defining lambda

Mrk_sub <- Mrk_cln[which(rownames(Mrk_cln) %in% as.character(pheno_dat$Pedigree)), ] - 1
pheno_sub <- pheno_dat %>% 
  filter(as.character(Pedigree) %in% rownames(Mrk_sub))

nIter <- 20000
burnIn <- 2000

# Bayesian LASSO
ETA <- list(MRK = list(X = Mrk_sub, model = 'BL'))
fmBL <- BGLR(y = pheno_sub$dEBV, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = 'Pepper/BL_')

mu <- fmBL$ETA$MRK$b

D_BL <- diag(1 / ncol(Mrk_cln) * (2 * allele_p * (1 - allele_p) * mu * mu))
colnames(D_BL) <- colnames(Mrk_sub)
rownames(D_BL) <- colnames(Mrk_sub)

G_BL <- Z %*% D_BL %*% t(Z)

G_BL_sub <- G_BL[which(rownames(G_BL) %in% as.character(pheno_dat$Pedigree)),
                 which(colnames(G_BL) %in% as.character(pheno_dat$Pedigree))]


# Question, how to define lambda and check convergence of Bayesian model?

# G inverse for Asreml with synbreed (write.realtionshipMatrix)

res_1 <- write.relationshipMatrix(G_BL_sub, sorting = 'ASReml',
                                type = c('ginv'), digits = 10)
names(res_1) <- c('row', 'column', 'coefficient')


# res <- KinshipMatrix(G_BL_sub)

# Test G inverse in asreml

# Bayesian LASSO with user defined prior

# Prior hyperparameter values
# sigmaE2, residual variance
r <- 1-0.36 # 0.36 is heritability of average fruit wall thickness
mode.sigE <- r*var(pheno_sub$dEBV)
dfe <- 3
Se <- mode.sigE*(dfe + 2)

# lambda
mode.sigL <- (1 - r)*var(pheno_sub$dEBV)
rate <- 2*Se/mode.sigL*sum(colMeans(Mrk_sub)^2)
shape <- 1.1

# Set priors
prior <- list(varE = list(S0 = Se, df0 = dfe),
              lambda = list(type = 'random',
                            shape = shape, 
                            rate = rate))
# lambda <- (1-0.3)/0.3 # heritability is 0.3
# 
# prior <- list(lambda = list(type = 'fixed', 
#                             value = lambda))

ETA <- list(MRK = list(X = Mrk_sub, model = 'BL', prior))
fmR_l <- BGLR(y = pheno_sub$dEBV, 
              ETA = ETA, 
              nIter = nIter, 
              burnIn = burnIn, 
              saveAt = 'Pepper/BL_lambda__')

mu <- fmR_l$ETA$MRK$b

D_BL <- diag(1 / ncol(Mrk_cln) * (2 * allele_p * (1 - allele_p) * mu * mu))
colnames(D_BL) <- colnames(Mrk_sub)
rownames(D_BL) <- colnames(Mrk_sub)

G_BL <- Z %*% D_BL %*% t(Z)

# G_BL_sub <- G_BL[which(rownames(G_BL) %in% pheno_dat$Pedigree), which(colnames(G_BL) %in% pheno_dat$Pedigree)]

G_BL_sub <- G_BL[which(rownames(G_BL) %in% as.character(pheno_sub$Pedigree)), 
                 which(colnames(G_BL) %in% as.character(pheno_sub$Pedigree))]

res_1 <- write.relationshipMatrix(G_BL_sub, sorting = 'ASReml',
                                  type = c('ginv'), digits = 10)
names(res_1) <- c('row', 'column', 'coefficient')

# visualize G matrix
library(ggcorrplot)
ggcorrplot(G_base_sub[,ncol(G_BL_sub):1]) # inverse columns to get  change the orientation of diagonal of G matrix

############################ GLBUP with Asreml #####################

library(asreml)

gblup_mdl <- asreml(fixed = dEBV ~ 1, 
                    random = ~giv(Pedigree, var = T),
                    ginverse = list(Pedigree = res_1),
                    rcov = ~units,
                    data = pheno_sub,
                    na.method.X = 'include')
pred_gblup <- predict(gblup_mdl, classify = 'Pedigree')$predictions$pval

cor(pred_gblup$predicted.value, pheno_sub$dEBV)
plot(pred_gblup$predicted.value, pheno_sub$dEBV)

