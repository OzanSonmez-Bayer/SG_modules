# Elastic Net

library(aws.s3)
library(tidyverse)
library(caret)
library(glmnet)
library(jsonlite)
library(data.table)
library(Rcpp)


### Elastic net with random search for alpha 

TrainModels_randomsearch <- function(pheno, geno, trt, outdir){
  pheno_geno <- pheno %>% 
    dplyr::rename(Pedigree =  PEDIGREE_NAME) %>% 
    filter(trait == trt) %>% 
    inner_join(geno,  by = 'Pedigree')
  
  # partition data
  
  set.seed(2020)
  
  training_indx <- createDataPartition(pheno_geno$dEBV, times = 1, p = 0.9, list = FALSE)
  
  training_set <- pheno_geno[training_indx, ]
  testing_set <- pheno_geno[-training_indx, ]
  
  y_train <- training_set %>% 
    select(dEBV) %>% 
    as.matrix()
  X_train <- training_set %>% 
    select(-c(dEBV, Pedigree, trait, ProgenyGermID)) %>% 
    as.matrix()
  
  y_test <- testing_set %>% 
    select(dEBV) %>% 
    as.matrix()
  X_test <- testing_set %>% 
    select(-c(dEBV, Pedigree, trait, ProgenyGermID)) %>% 
    as.matrix()
  
  # Cross validation with training set to tune hyper-parameters
  
  train_control <- trainControl(method = "repeatedcv",
                                number = 5,
                                repeats = 1,
                                search = "random",
                                verboseIter = TRUE)
  
  # Train the model: Elastic net
  elastic_net_model <- train(dEBV ~ .,
                             data = cbind(y_train, X_train),
                             method = "glmnet",
                             tuneLength = 50,
                             trControl = train_control)
  
  elastic_net_model$bestTune
  
  
  # LASSO 
  # lambda <- 10^seq(-5, 5, length = 100)
  # lasso_model <- train(
  #   dEBV ~ .,
  #   data = cbind(y_train, X_train),
  #   method = "glmnet",
  #   trControl = trainControl("cv", number = 5),
  #   tuneGrid = expand.grid(alpha = 1, lambda = lambda)
  # )
  
  lasso_tuning <- cv.glmnet(X_train, y_train, nfolds = 5, alpha = 1)
  lasso_model <- glmnet(X_train, y_train, lambda = lasso_tuning$lambda.min, alpha = 1)
  
  # Ridge Regression
  # 
  # ridge_model <- train(
  #   dEBV ~ .,
  #   data = cbind(y_train, X_train),
  #   method = "glmnet",
  #   trControl = trainControl("cv", number = 5),
  #   tuneGrid = expand.grid(alpha = 0, lambda = lambda)
  # )
  
  ridge_tuning <- cv.glmnet(X_train, y_train, nfolds = 5, alpha = 0)
  ridge_model <- glmnet(X_train, y_train, lambda = lasso_tuning$lambda.min, alpha = 0)
  
  # Save finalized model
  saveRDS(list(enet = elastic_net_model, lasso = lasso_model, ridge = ridge_model), 
          paste0(outdir, trt, '.rds'))
  
  accuracy_df <- evaluationmetric(y_test, X_test, trt, elastic_net_model, lasso_model, ridge_model)
  
  return(accuracy_df)
}

# Elastic Net with pre-defined alpha 

TrainModels_predefined <- function(pheno, geno, trt, alpha_list, outdir){
  
  set.seed(20200416)
  
  results_df <- data.frame()
  marker_list <- list()
  
  pheno_geno <- pheno %>% 
    dplyr::rename(Pedigree =  PEDIGREE_NAME) %>% 
    filter(trait == trt) %>% 
    inner_join(geno,  by = 'Pedigree')
  
  training_indx <- createDataPartition(pheno_geno$dEBV, times = 1, p = 0.9, list = FALSE)
  
  training_set <- pheno_geno[training_indx, ]
  testing_set <- pheno_geno[-training_indx, ]
  
  y_train <- training_set %>% 
    select(dEBV) %>% 
    as.matrix()
  X_train <- training_set %>% 
    select(-c(dEBV, Pedigree, trait, ProgenyGermID)) %>% 
    as.matrix()
  
  y_test <- testing_set %>% 
    select(dEBV) %>% 
    as.matrix()
  X_test <- testing_set %>% 
    select(-c(dEBV, Pedigree, trait, ProgenyGermID)) %>% 
    as.matrix()
  
  
  for(par_i in 1:length(alpha_list)){
    mdl_tuning <- cv.glmnet(X_train, y_train, nfolds = 5, alpha = alpha_list[par_i])
    mdl <- glmnet(X_train, y_train, lambda = mdl_tuning$lambda.min, alpha = alpha_list[par_i])
    
    y_hat <- mdl %>% predict(X_test) %>% as.vector()
    rsq <- cor(y_test, y_hat)^2
    
    mdl_result <- data.frame('Alpha' = alpha_list[par_i], 
                             'Accuracy' = round(rsq, 3),
                             'Number of Marker' = length(coef(mdl)@i) - 1)
    results_df <- rbind(results_df, mdl_result)
    marker_list <- c(marker_list, assign(paste0('alpha_', alpha_list[par_i]), list(coef(mdl)@Dimnames[[1]][coef(mdl)@i[-1]])))
    
  }
  
  saveRDS(list(accuracy = results_df, marker = marker_list), 
          paste0(outdir, trt, '_ENET.rds'))
  
}


evaluationmetric <- function(y, x, trtname, enet, lasso, ridge){
  y_hat_enet <- predict(enet, x)
  rsq_enet <- cor(y, y_hat_enet)^2
  
  y_hat_lasso <- lasso %>% 
    predict(x) %>% 
    as.vector()
  rsq_lasso <- cor(y, y_hat_lasso)^2
  
  y_hat_ridge <- ridge %>% 
    predict(x) %>% 
    as.vector()
  rsq_ridge <- cor(y, y_hat_ridge)^2
  
  return(data.frame('Trait' = trtname,
                    'ENET' = round(rsq_enet, 3),
                    'ENET_marker' = length(coef(enet$finalModel, enet$finalModel$lambdaOpt)@i) - 1,
                    'LASSO' = round(rsq_lasso, 3),
                    'LASSO_marker' = length(coef(lasso)@i) - 1,
                    'Ridge'= round(rsq_ridge, 3), 
                    'Ridge_marker' = length(coef(ridge)@i) - 1))
}


########## Train Elastic Model for hybrids ###########

trainingModel <- function(pheno_raw, geno, trait, alpha_list, outputfile, validation_only, testing_season){
  
  # Training data
  pheno <- pheno_raw %>% 
    filter(OBSRVTN_REF_CD == trait & GROWSEASON != testing_season) %>% 
    select(PEDIGREE_NAME, TRAIT_VALUE) %>% 
    mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
    group_by(PEDIGREE_NAME) %>% 
    summarize(lsmean = mean(TRAIT_VALUE, na.rm = TRUE))
  
  head(pheno)
  
  pheno_geno <- pheno %>% 
    inner_join(geno, by = c('PEDIGREE_NAME' = 'Pedigree'))
  
  set.seed(2020)
  
  training_indx <- createDataPartition(pheno_geno$lsmean, times = 1, p = 0.9, list = FALSE)
  
  training_set <- pheno_geno[training_indx, ]
  testing_set <- pheno_geno[-training_indx, ]
  
  y_train <- training_set %>% 
    select(lsmean) %>% 
    as.matrix()
  X_train <- training_set %>% 
    select(-c(lsmean, PEDIGREE_NAME, GermID)) %>% 
    as.matrix()
  
  y_test <- testing_set %>% 
    select(lsmean) %>% 
    as.matrix()
  X_test <- testing_set %>% 
    select(-c(lsmean, PEDIGREE_NAME, GermID)) %>% 
    as.matrix()
  
  # True testing data
  if (validation_only == 'True'){
    pheno_testing <- pheno_raw %>% 
      filter(OBSRVTN_REF_CD == trait & GROWSEASON == testing_season) %>% 
      select(PEDIGREE_NAME, TRAIT_VALUE) %>% 
      mutate(TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
      group_by(PEDIGREE_NAME) %>% 
      summarize(lsmean = mean(TRAIT_VALUE, na.rm = TRUE))
    
    pheno_geno_testing <- pheno_testing %>% 
      inner_join(geno, by = c('PEDIGREE_NAME' = 'Pedigree'))
    
    y_true_testing <- pheno_geno_testing %>% 
      select(lsmean) %>% 
      as.matrix()
    
    X_true_testing <- pheno_geno_testing %>% 
      select(-c(lsmean, PEDIGREE_NAME, GermID)) %>% 
      as.matrix()
  }else{
    print('Only validation')
  }
 
  
  rm(pheno)
  rm(geno)
  gc()
  
  
  results_df <- data.frame()
  for(par_i in 1:length(alpha_list)){
    print(alpha_list[par_i])
    mdl_tuning <- cv.glmnet(X_train, y_train, nfolds = 5, alpha = alpha_list[par_i])
    mdl <- glmnet(X_train, y_train, lambda = mdl_tuning$lambda.min, alpha = alpha_list[par_i])
    
    y_hat <- mdl %>% predict(X_test) %>% as.vector()
    rsq <- cor(y_test, y_hat)
    
    if(validation_only == 'True'){
      y_hat_testing <- mdl %>% predict(X_true_testing) %>% as.vector()
      rsq_testing <- cor(y_true_testing, y_hat_testing)
      
      mdl_result <- data.frame('Alpha' = alpha_list[par_i], 
                               'Validation Accuracy' = round(rsq, 3),
                               'Testing Accuracy' = round(rsq_testing,3),
                               'Number of Marker' = length(coef(mdl)@i) - 1)
    }else{
      print('Only validation-2')
      mdl_result <- data.frame('Alpha' = alpha_list[par_i], 
                               'Validation Accuracy' = round(rsq, 3),
                               'Number of Marker' = length(coef(mdl)@i) - 1)
    }
    
    results_df <- rbind(results_df, mdl_result)
    # marker_list <- c(marker_list, assign(paste0('alpha_', alpha_list[par_i]), list(coef(mdl)@Dimnames[[1]][coef(mdl)@i[-1]])))
    
  }
  saveRDS(results_df, file = paste0(outputfile, trait, '.rds'))
  
}
