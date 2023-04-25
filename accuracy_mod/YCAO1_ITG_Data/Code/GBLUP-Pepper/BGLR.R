set.seed(20200217)

library(BGLR)
library(coda)
# require(devtools)
# install_version("MCMCpack", version = "1.4-0", repos = "http://cran.us.r-project.org")

library(MCMCpack)

# BGLR requires matrix of markers
# marker matrix is centered on zero for ease of calculation

geno_centered <- Mrk_sub_1 - 1

y <- as.array(pred_df$dEBV)

samplesize <- ceiling(53 * 0.2)
fold <- 2

whichNa <- sample(1:length(y), size = samplesize, replace = FALSE)
yNa <- y
yNa[whichNa] <- NA

# Basic Bayesian LASSO
fmR <- BGLR(y = yNa, response_type = 'gaussian', a = NULL, b = NULL, 
            ETA =  list(MRK = list(X = geno_centered, model = 'BL')),
            nIter = 2000, 
            burnIn = 100, 
            thin  =10)

summary(fmR)

# --------------------> Summary of data & model <-------------------- 
#   
# Number of phenotypes= 42 
# Min (TRN)=  -1.445226 
# Max (TRN)=  1.758436 
# Variance of phenotypes (TRN)= 0.5422  # heritability  = 0.5422/(0.5422 + 0.2943)
# Residual variance= 0.2943 
# N-TRN= 42  N-TST= 11 
# Correlation TRN= 0.9747 
# 
# -- Linear Predictor -- 
#   
#   Intercept included by default
# Coefficientes in ETA[ 1 ] ( MRK ) modeled as in  BL 
# 
# ------------------------------------------------------------------


# Calculate model statistics
df1 <- data.frame(fmR$yHat[whichNa], y[whichNa])
PredAbi <- cor(df1[,2], df1[,1])
Rank <- cor(df1[,2], df1[,1], method = 'spearman')
MSE <- mean((df1[,2] - df1[, 1])^2)
Bias <- coef(lm(y[whichNa] ~ fmR$yHat[whichNa]))[2]
table_accuracy <- cbind(PredAbi, Rank, MSE, Bias)
table_accuracy


# Other Bayesian Methods

nIter <- 20000
burnIn <- 2000

# Bayesian LASSO
ETA <- list(MRK = list(X = geno_centered, model = 'BL'))
fmBL <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = 'BL_')

## Bayesian Ridge Regression (Gaussian Prior)
ETA <- list(MRK = list(X = geno_centered, model = 'BRR'))
fmBRR <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, 
              saveAt = 'BRR_')


## Bayes A (Scaled-t prior)
ETA$MRK$model <- 'BayesA'
fmBA <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, 
              saveAt = 'BA_')

## Bayes B (point of mass at zero + scaled-t slab)
ETA$MRK$model <- 'BayesB'
fmBB <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, 
             saveAt = 'BB_')


## Bayes C (point of mass at zero + scaled-t slab)
ETA$MRK$model <- 'BayesC'
fmBC <- BGLR(y = yNa, ETA = ETA, nIter = nIter, burnIn = burnIn, 
             saveAt = 'BC_')



### Calcualte correlations
r_BL <- cor(y[whichNa], fmBL$yHat[whichNa])
r_BRR <- cor(y[whichNa], fmBRR$yHat[whichNa])
r_BA <- cor(y[whichNa], fmBA$yHat[whichNa])
r_BB <- cor(y[whichNa], fmBB$yHat[whichNa])
r_BC <- cor(y[whichNa], fmBC$yHat[whichNa])

table_accuracy <- data.frame(rbind(r_BL, r_BRR, r_BA, r_BB, r_BC),
                             rbind(fmBL$varE, fmBRR$varE, fmBA$varE, fmBB$varE, fmBC$varE),
                             rbind(fmBL$fit$pD, fmBRR$fit$pD, fmBA$fit$pD, fmBB$fit$pD, fmBC$fit$pD), 
                             rbind(fmBL$fit$DIC, fmBRR$fit$DIC, fmBA$fit$DIC, fmBB$fit$DIC, fmBC$fit$DIC))

colnames(table_accuracy) <- c('PredAbi', 'varE', 'pD', 'DIC')
rownames(table_accuracy) <- c('BLASSO', 'BRidge', 'BayersA', 'BayersB', 'BayersC')

round(table_accuracy, 3)
