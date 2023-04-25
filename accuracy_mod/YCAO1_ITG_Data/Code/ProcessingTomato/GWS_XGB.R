# GBLUP for Processing Tomato

library(readxl)
library(tidyverse)
library(ggplot2)
library(aws.s3)
library(jsonlite)
library(Matrix)
# install.packages('xgboost')
library(xgboost)
library(caret)


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


procTom_taqman_2[1:5,1:5]

taqman_full <- procTom_taqman_2 %>% 
  dplyr::select(-Inventory_SampleID) %>% 
  unique()
rownames(taqman_full) <- taqman_full$Pedigree
taqman_full <- taqman_full %>% 
  dplyr::select(-Pedigree)

taqman_full_1 <- gen2Additive(taqman_full)
rownames(taqman_full_1) <- rownames(taqman_full)

G_rrBLUP <- rrBLUP::A.mat((taqman_full_1 - 1), 
                          min.MAF = 0.01,  
                          impute.method = 'mean', 
                          tol = 0.02, 
                          return.imputed = TRUE)

G_imputed <- G_rrBLUP$imputed
hist(diag(G_rrBLUP$A))
range(diag(G_rrBLUP$A)) 



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

# MAT, FRFRMH, PLTVG, PLCMP, FQUAL, PYLDPA, SAMPUW, OST, AVJB, LBRIX, FBRIX, AVGHBRIX, FZUNI, MATR, EARLY

hyperparam_list <- list()

trt <- 'PLTVG'
pheno_df <- pheno %>% 
  filter(trait == trt) %>% 
  mutate(dEBV = predicted.value/reliability)

length(intersect(as.character(pheno_df$PEDIGREE_NAME), rownames(G_imputed))) # 496 inbreds have been genotyped

pheno_sub <- pheno_df %>% 
  filter(as.character(PEDIGREE_NAME) %in% rownames(G_imputed)) %>%
  dplyr::select(PEDIGREE_NAME, dEBV)
pheno_sub <- pheno_sub[!duplicated(pheno_sub$PEDIGREE_NAME),]
rownames(pheno_sub) <- pheno_sub$PEDIGREE_NAME
pheno_sub$PEDIGREE_NAME <- as.character(pheno_sub$PEDIGREE_NAME)
# pheno_sub <- pheno_sub %>% 
#   dplyr::select(-PEDIGREE_NAME)


G_sub <- G_imputed[rownames(G_imputed) %in% pheno_sub$PEDIGREE_NAME,]
G_sub <- G_sub[pheno_sub$PEDIGREE_NAME,]

######################## Training - Tessting Split ##################
set.seed(20200422)
trainIndex <- createDataPartition(pheno_sub$dEBV, p = 0.8, list = FALSE)
x_train <- G_sub[trainIndex, ]
x_test <- G_sub[-trainIndex, ]

y_train <- pheno_sub$dEBV[trainIndex]
y_test <- pheno_sub$dEBV[-trainIndex]


###################### Automatic Parameter tuning with Caret (Random Search) ##############
dt <- data.frame(dEBV = y_train, x_train)

fitControl <- trainControl(method = "cv",
                           number = 3,
                           search = "random")

mdl_xgb <- train(dEBV ~. , data = dt, 
                 method = "xgbTree", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 eval.metric = "rmse",
                 verbose = FALSE)

cor(y_test, predict(mdl_xgb, newdata = x_test))
cor(y_train, predict(mdl_xgb, nenewdata = x_train))

hyperparam_list[[trt]] <- mdl_xgb$bestTune

######################## Fine-Tune Hyperparameters with Grid Search  ##################

# nrounds
# eta: learning rate
# max_depth
# gamma
# colsample_bytree
# subsample

nrounds <- 200

# step 1: number of iterations and the learning rate
tune_grid <- expand.grid(
  nrounds = seq(from = 1, to = nrounds, by = 10),
  eta = c(0.025, 0.05, 0.1, 0.3),
  max_depth = c(2, 3, 4, 5, 6),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)


tune_control <- caret::trainControl(
  method = 'cv',
  number = 3,
  verboseIter = FALSE, 
  allowParallel = TRUE
)


xgb_tune <- caret::train(
  x = x_train, 
  y = y_train,
  trControl = tune_control, 
  tuenGrid = tune_grid,
  method = 'xgbTree',
  verbose = TRUE
)


tuneplot <- function(x, probs = .90) {
  ggplot(x) +
    coord_cartesian(ylim = c(quantile(x$results$RMSE, probs = probs), min(x$results$RMSE))) +
    theme_bw()
}
# 
# 
tuneplot(xgb_tune)

xgb_tune$bestTune
# nrounds max_depth eta gamma colsample_bytree min_child_weight subsample
#  50         3     0.3     0              0.6                1         1
min(xgb_tune$results$RMSE) # 30.0151

# ggplot(xgb_tune)


# step 2: Max depth and min child weight
tune_grid2 <- expand.grid(
  nrounds = seq(from = 10, to = nrounds, by = 5),
  eta = xgb_tune$bestTune$eta,
  max_depth = ifelse(xgb_tune$bestTune$max_depth == 2,
                     c(xgb_tune$bestTune$max_depth:4),
                     xgb_tune$bestTune$max_depth - 1:xgb_tune$bestTune$max_depth + 1),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(1, 2, 3),
  subsample = 1
)

xgb_tune2 <- caret::train(
  x = x_train,
  y = y_train,
  trControl = tune_control,
  tuneGrid = tune_grid2,
  method = "xgbTree",
  verbose = TRUE
)


tuneplot(xgb_tune2)

xgb_tune2$bestTune
# nrounds max_depth eta gamma colsample_bytree min_child_weight subsample
# 20         3      0.3     0                1                3         1
min(xgb_tune2$results$RMSE) # 30.46265


# Step 3: column and row sampling
tune_grid3 <- expand.grid(
  nrounds = seq(from = 10, to = nrounds, by = 5),
  eta = xgb_tune$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = 0,
  colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = c(0.5, 0.75, 1.0)
)

xgb_tune3 <- caret::train(
  x = x_train,
  y = y_train,
  trControl = tune_control,
  tuneGrid = tune_grid3,
  method = "xgbTree",
  verbose = TRUE
)

tuneplot(xgb_tune3, probs = .95)

xgb_tune3$bestTune
# nrounds max_depth eta gamma colsample_bytree min_child_weight subsample
#   30         3    0.3     0              0.4                3         1

min(xgb_tune3$results$RMSE) # 29.845

# step 4: Gamma

tune_grid4 <- expand.grid(
  nrounds = seq(from = 10, to = nrounds, by = 5),
  eta = xgb_tune$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0),
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune4 <- caret::train(
  x = x_train,
  y = y_train,
  trControl = tune_control,
  tuneGrid = tune_grid4,
  method = "xgbTree",
  verbose = TRUE
)

tuneplot(xgb_tune4)

xgb_tune4$bestTune

# nrounds max_depth eta gamma colsample_bytree min_child_weight subsample
#  20         3     0.3   0.5              0.4                3         1

min(xgb_tune4$results$RMSE) # 30.102

# Step 5: Reduce Learning Rate 
tune_grid5 <- expand.grid(
  nrounds = seq(from = 10, to = 100, by = 5),
  eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = xgb_tune4$bestTune$gamma,
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune5 <- caret::train(
  x = x_train,
  y = y_train,
  trControl = tune_control,
  tuneGrid = tune_grid5,
  method = "xgbTree",
  verbose = TRUE
)

tuneplot(xgb_tune5)

xgb_tune5$bestTune
min(xgb_tune5$results$RMSE) # 29.395

# Fitting Model
(final_grid <- expand.grid(
  nrounds = xgb_tune5$bestTune$nrounds,
  eta = xgb_tune5$bestTune$eta,
  max_depth = xgb_tune5$bestTune$max_depth,
  gamma = xgb_tune5$bestTune$gamma,
  colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
  min_child_weight = xgb_tune5$bestTune$min_child_weight,
  subsample = xgb_tune5$bestTune$subsample
))

# nrounds eta max_depth gamma colsample_bytree min_child_weight subsample
#  95     0.1         3   0.5              0.4                3         1

x_cbi <- rbind(x_train, x_test)
y_cbi <- c(y_train, y_test)

train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

(xgb_model <- caret::train(
  x = x_train,
  y = y_train,
  trControl = train_control,
  tuneGrid = final_grid,
  method = "xgbTree",
  verbose = TRUE
))


(xgb_model_rmse <- ModelMetrics::rmse(y_test, predict(xgb_model, newdata = x_test)))
(xgb_model_rmse <- cor(y_test, predict(xgb_model, newdata = x_test)))
cor(y_train, predict(xgb_model, newdata = x_train))


hyperparam_list <- list()
hyperparam_list[['MAT']] <- xgb_model$bestTune


################################# XGBOOST #########################


dtrain <- xgb.DMatrix(data = x_train, label = y_train)
dtest <- xgb.DMatrix(data = x_test, label = y_test)
# watchlist <- list(train = dat_1, test = dtest)

param <- list(max_depth = 3, 
             eta = 0.1, 
             nthread = 2, 
             gamma = 0.5,
             colsample_tytree = 0.4, 
             min_child_weight = 3, 
             subsample = 1)
mdl_xgb <- xgb.train(param, 
                     dtrain, 
                     eval.metric = "rmse",
                     nrounds = 95)

pred_xgb <- predict(mdl_xgb, x_test)
cor(y_test, pred_xgb)

cor(y_train, predict(mdl_xgb, x_train))

xgb_vimp <- xgb.importance(model = mdl_xgb) 
xgb.plot.importance(xgb_vimp)

