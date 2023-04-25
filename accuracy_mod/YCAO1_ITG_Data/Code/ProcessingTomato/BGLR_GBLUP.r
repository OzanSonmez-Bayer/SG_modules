# BGLR
library(aws.s3)
library(httr)
install.packages('BGLR')
library(BGLR)
library(tidyverse)

source('vaultCredentials.R')

VaultToken <- content(POST(url = paste0(VaultURL, "auth/approle/login"), body= AppRoleCredsJson, encode = 'json'))
print(VaultToken)
AWSCreds <- content(GET(url = paste0(VaultURL, VaultSecretPath), httr::add_headers('X-Vault-Token' = VaultToken$auth$client_token)))
print(AWSCreds)
Sys.setenv("AWS_ACCESS_KEY_ID" = AWSCreds$data$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = AWSCreds$data$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION" = "us-east-1")

geno <- s3readRDS(object = 'Tomato/ProcessingTomato/9Z/geno_all.rds', 
                  bucket = 'veg-apd-sdi-predictiveanalytics-prod-geno-data')

pheno <- s3readRDS(object = 'ycao1/ITG/Tomato/9Z/pheno_training.rds',
                   bucket = 'veg-apd-sdi-predictiveanalytcs-prod-workspace')

### Center G matrix 

geno_centered <- geno[, 3:ncol(geno)]
rownames(geno_centered) <- geno$Pedigree
geno_centered$pedigree <- rownames(geno_centered)

traitlist <- c('AVJB', 'FQUAL', 'FRFRMH', 'FZUNI', 'LBRIX', 'MAT', 'MATR', 'OST', 'PLCMP', 'PLTVG', 'PYLDPA', 'SAMPUW')

model_results <- data.frame(Trait = rep(traitlist,rep(5,length(traitlist))), 
                            Model = rep(c('Bayesian LASSO', 'Bayesian Ridge', 'BayesA', 'BayesB', 'BayesC'), length(traitlist)),
                            Accuracy = NA, 
                            pD = NA, 
                            DIC = NA) 

traitlist <- c('MAT', 'PLCMP', 'PLTVG')
  
for (trait in traitlist){
  
  print(trait)
  # trait <- 'AVGHBRIX'
  
  pheno_sub <- pheno %>% 
    filter(OBSRVTN_REF_CD == trait) %>% 
    rename(pedigree = PEDIGREE_NAME) %>% 
    mutate(TREP = paste(REP_NUMBER, TRACK_ID,sep='_'),
           TRAIT_VALUE = as.numeric(TRAIT_VALUE)) %>% 
    group_by(pedigree) %>% 
    summarise(lsm = mean(TRAIT_VALUE))
  
  pheno_geno <- pheno_sub %>% 
    inner_join(geno_centered, by = 'pedigree')
  
  # mask testing season
  yNA <- pheno_geno$lsm
  testing_ped <- unique(pheno$PEDIGREE_NAME[which(pheno$GROWSEASON == '2019:04' & pheno$OBSRVTN_REF_CD == trait)])
  yNA[which(pheno_geno$pedigree %in% testing_ped)] <- NA
  
  ETA <- list(MRK = list(X = pheno_geno[,3:nrow(pheno_geno)], model = 'BL'))
  
  ptm <- proc.time()
  fm <- BGLR(y = yNA, ETA = ETA, nIter = 12000, burnIn = 2000, saveAt = '/mnt/GWS/ProcTomato/BGLR/' )
  proc.time() - ptm
  
  save(fm, file = paste0('/mnt/GWS/ProcTomato/BGLR/', trait, '_BL.rds'))
  
  # Predictions 
  yHat <- fm$yHat
  # temp <- range(c(pheno_geno$lsm, yHat))
  # plot(yHat ~ pheno_geno$lsm, xlab = 'Observed', ylab = 'Predicted', col = 2, 
       # xlim = temp, ylim = temp)
  # testing accuracy 
  r_BL <- cor(pheno_geno$lsm[which(is.na(yNA))], yHat[which(is.na(yNA))], use = 'complete.obs')
  
  # Bayesian Ridge
  ETA$MRK$model <- 'BRR'
  fmBRR <- BGLR(y = yNA, ETA = ETA, nIter = 12000, burnIn = 2000, saveAt = '/mnt/GWS/ProcTomato/BGLR/' )
  save(fmBRR, file =  paste0('/mnt/GWS/ProcTomato/BGLR/', trait, '_', ETA$MRK$model,'.rds'))
  
  yHat_BRR <- fmBRR$yHat
  # temp <- range(c(pheno_geno$lsm, yHat_BRR))
  # plot(yHat_BRR ~ pheno_geno$lsm, xlab = 'Observed', ylab = 'Predicted', col = 2, 
  #      xlim = temp, ylim = temp)
  # testing accuracy 
  r_BRR <- cor(pheno_geno$lsm[which(is.na(yNA))], yHat_BRR[which(is.na(yNA))], use = 'complete.obs')
  
  
  # BayesA
  
  ETA$MRK$model <- 'BayesA'
  fmBA <- BGLR(y = yNA, ETA = ETA, nIter = 12000, burnIn = 2000, saveAt = '/mnt/GWS/ProcTomato/BGLR/' )
  save(fmBA, file =  paste0('/mnt/GWS/ProcTomato/BGLR/', trait, '_', ETA$MRK$model,'.rds'))
  
  yHat_BA <- fmBA$yHat
  r_BA <- cor(pheno_geno$lsm[which(is.na(yNA))], yHat_BA[which(is.na(yNA))], use = 'complete.obs')
  
  # BayesB
  
  ETA$MRK$model <- 'BayesB'
  fmBB <- BGLR(y = yNA, ETA = ETA, nIter = 12000, burnIn = 2000, saveAt = '/mnt/GWS/ProcTomato/BGLR/' )
  save(fmBB, file =  paste0('/mnt/GWS/ProcTomato/BGLR/', trait, '_', ETA$MRK$model,'.rds'))
  
  yHat_BB <- fmBB$yHat
  r_BB <- cor(pheno_geno$lsm[which(is.na(yNA))], yHat_BB[which(is.na(yNA))], use = 'complete.obs')
  
  # BayesC
  
  ETA$MRK$model <- 'BayesC'
  fmBC <- BGLR(y = yNA, ETA = ETA, nIter = 12000, burnIn = 2000, saveAt = '/mnt/GWS/ProcTomato/BGLR/' )
  save(fmBC, file =  paste0('/mnt/GWS/ProcTomato/BGLR/', trait, '_', ETA$MRK$model,'.rds'))
  
  yHat_BC <- fmBC$yHat
  r_BC <- cor(pheno_geno$lsm[which(is.na(yNA))], yHat_BC[which(is.na(yNA))], use = 'complete.obs')
  
  temp_acc <- model_results[which(model_results$Trait == trait),]
  temp_acc$Accuracy <- c(r_BL, r_BRR, r_BA, r_BB, r_BC)
  temp_acc$pD <- c(fm$fit$pD, fmBRR$fit$pD, fmBA$fit$pD, fmBB$fit$pD, fmBC$fit$pD)
  temp_acc$DIC <- c(fm$fit$DIC, fmBRR$fit$DIC, fmBA$fit$DIC, fmBB$fit$DIC, fmBC$fit$DIC)
  model_results[which(model_results$Trait == trait), ] <- temp_acc
}


# Metrics to compare models 
# fm$fit$pD # the effetive number of predictors, shwoing the complexity of the model
# fm$fit$DIC # deviance inoformation criteria. Smaller DIC values are preferred. 

# Visualize accuracy, pD, and DIC

g1 <- ggplot(model_results, aes(x=Model, y=Accuracy, group=Trait)) +
  geom_line(aes(color=Trait))+
  geom_point(aes(color=Trait)) + 
  ggtitle('Accuracy') + 
  theme_classic() +
  theme(legend.position = "none")


g2 <- ggplot(model_results, aes(x=Model, y=pD, group=Trait)) +
  geom_line(aes(color=Trait))+
  geom_point(aes(color=Trait)) + 
  ggtitle('Model Complexity') + 
  theme_classic() +
  theme(legend.position = "none")

g3 <- ggplot(model_results, aes(x=Model, y=DIC, group=Trait)) +
  geom_line(aes(color=Trait))+
  geom_point(aes(color=Trait)) + 
  ggtitle('DIC: Deviance Information Criteria') + 
  theme_classic() 

ggarrange(g1, g2, g3, 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

ggarrange(                                              
          ggarrange(g1, g2, ncol = 2, labels = c("", "")), 
          g3,
          nrow = 2, 
          labels = ""                                        
)
