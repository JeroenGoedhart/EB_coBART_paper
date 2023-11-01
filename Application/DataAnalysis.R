setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Application")


### loading relevant libraries ###
library(dbarts)
library(loo)
library(pROC)

gc()


load("DatForPaper.Rdata")
Xtrain = dat$Xtrain # feature matrix, 140 features of which 67 copy number variation, 69 mutation, 3 translocation, and the IPI
Ytrain = dat$Ytrain # two year progression free survival (yes = 0, no = 1)
CoData = dat$CoData # Codata matrix, consists of p values (on logit scale) of each feature (Estimated on independent cohort with same outcome), and feature type
Xtest = dat$Xtest
Ytest = dat$Ytest


remove(dat); gc()




###############################################################################
################################ Fit EB-coBART ################################
##############################################################################
source('sidefunctions.R') # load side functions consists of
# EstimateLeafNode: empirical bayes estimates of k hyperparameter (leaf node) of BART
# FiniteSum: only sum observations that have finite values
# getDepth: computes the depth of each node of all trees of a BART fit
# LikelihoodBin: computes likelihood values for binary variables based on bernoulli distribution
# LikelihoodTreestructure: compute the prior probability on log scale of a tree structure to occur as a function of
                            # alpha and beta hyperparameters
getwd()
## Initialization parameters
nIter = 18 # number of iterations that EBcoBART runs
probs <- rep(1/ncol(Xtrain),ncol(Xtrain)) # initial covariate weights, each variable gets the same  weight of 1/p


# initialize model in rigid or flexible tree setting
#model <- "rigid"
model <- "flexible"

EBupdates <- T #set to false if not updating alpha and k

if (model=="rigid"){k = 1; base = .1; power = 4.0}
if (model=="flexible"){k=2; base = .95; power = 2.0}


# storage containers
EstimatedProbs <- matrix(NA, nrow = nIter+1, ncol = ncol(Xtrain))
Codatamodels <- vector("list", length = nIter)
EstimatedProbs[1,] <- probs
row.names(EstimatedProbs) <- c(paste("Iter",0:(nIter), sep = " "))

WAICVector <- vector("list", length = nIter)

if (EBupdates == T){
  k_Update <- c()
  base_update <- c()
}

gc()


for (i in 1:nIter) {
  print(paste("iteration",i,sep = " "))
  
  ## STEP 1: Fit BART model ##
  fit <- bart(x.train = Xtrain, y.train = Ytrain,
              #x.test = Xtest,
              ndpost = 12000L,                   # number of posterior samples
              nskip = 12000L,                    # number of "warmup" samples to discard
              nchain = 10L,                      # number of independent, parallel chains
              ntree = 50L,                       # number of trees per chain
              verbose = F,
              usequants = F,
              k = k, base = base, power = power, # hyperparameters tree
              splitprobs = probs,                # prob that variable is chosen for split
              keeptrees = EBupdates,             #set to True if updating alpha and k and to False if not
              combinechains = T,
              seed = 4*i^2+202+3*i) 
  
  ## STEP 2: Extract information from BART fit ##
  VarsUsed <- colSums(fit$varcount)  # count of each variable occuring in the splitting rules
  VarsUsed <- VarsUsed/sum(VarsUsed) # normalize count of each variable to probabilities = pure EB updates of hyperparameter S
  Ypred = pnorm(fit$yhat.train)      # predictions for the training response
  
  
  ## STEP 3: Fit co-data model ##
  coDataModel <- glm(VarsUsed ~ Groups + p.values,
                     data=CoData,family=quasibinomial) # the model
  
  Codatamodels[[i]] <- coDataModel
  probs <- predict(coDataModel, type="response", newdata = CoData) # estimating the co-data moderated estimates of hyperparameter S
  probs[is.na(probs)] <- 0
  EstimatedProbs[i+1,] <- probs
  
  ## STEP 4: Estimate the WAIC for convergence ##
  Ypred[which(Ypred==0)] <- .0000000000000001
  Ypred[which(Ypred==1)] <- .9999999999999999
  LogLikMatrix = LikelihoodBin(Ypred = Ypred, Y = Ytrain)
  WAICVector[[i]] <- waic(LogLikMatrix)
  
  ## Optional step: update other hyperparameters (alpha and k) of BART using Empirical Bayes ##
  if (EBupdates == T) {
    trees <- extract(fit,"trees",chainNum = c(1:10), sampleNum=c(1:200,1000:1200,1800:2000)) # tree structures, for computation, we only select a part of the tree
    
    # Update leaf node parameter k
    k <- EstimateLeafNode(Trees = trees, ntree = 50)[2]
    k_Update[i] <- k
    
    # Update tree structure parameter alpha, we keep beta fixed
    trees <- trees[c("tree", "sample", "chain", "n","var","value" )]
    trees$depth <- unname(unlist(by(trees, trees[,c("tree", "sample", "chain")], getDepth)))
    base <- optim(base,LikelihoodTreeStructure,beta = power, Trees=trees, method = 'Brent', lower = .00001, upper = .9999999)$par
    base_update[i] <- base
    remove(trees)
  }
  remove(fit)
}

results = matrix(NA, nrow = nIter, ncol =3)
for (i in 1:nIter){
  results[i,] = WAICVector[[i]]$estimates[,1]}
waics = c(results[,3])
remove(results)
plot(waics)

minimum = function (x){
  for (i in 1:length(x)) {
    if (i==length(x)){
      return(i)
    }
    else {
      min = x[i]
      
      if (min<x[i+1]){
        return(i)
      }
    }
  }
}



min = minimum(waics)

probsmax <- EstimatedProbs[min,]
names(probsmax) <- names(Xtrain)
sum(probsmax)
probsmax <- probsmax/sum(probsmax)
sort(probsmax)
FinalCoDataModel <- Codatamodels[[min-1]]
summary(FinalCoDataModel)
k_min <- k_Update[min] 
base_min <- base_update[min]
hyperparameters = c(k_min, base_min, power)
names(hyperparameters) = c("k","base","power")

results = list(CoData = FinalCoDataModel, EstProbs = probsmax, WAICs = waics, hyperparameters = hyperparameters,
               alpha = base_update, k = k_Update)

getwd()
results = list(CoData = FinalCoDataModel, EstProbs = probsmax, WAICs = waics)

name <- paste("EBCoBARTfit",model,EBupdates,".Rdata", sep = "_")

save(results,file = name)

## load fitted EB-coBART models (TRUE: updates of alpha and k; FALSE: no updates of alpha and k)
load("Results/EBCoBARTfit_flexible_TRUE_.Rdata")
load("Results/EBCoBARTfit_flexible_FALSE_.Rdata")
load("Results/EBCoBARTfit_rigid_TRUE.Rdata")
load("Results/EBCoBARTfit_rigid_FALSE.Rdata")

load("DatForPaper.Rdata")
Xtrain = dat$Xtrain # feature matrix, 140 features of which 67 copy number variation, 69 mutation, 3 translocation, and the IPI
Ytrain = dat$Ytrain # two year progression free survival (yes = 0, no = 1)
Xtest = dat$Xtest
Ytest = dat$Ytest
remove(dat); gc()

probsmax = results$EstProbs # co-data moderated EB estimates of S
probsmax
Hyperparameters = results$hyperparameters # other hyperparmater estimates
k_min = Hyperparameters[1]; base_min = Hyperparameters[2]; power = Hyperparameters[3]
getwd()

library(dbarts)
library(pROC)
fit <- bart(x.train = Xtrain, y.train = Ytrain,
            x.test = Xtest,
            ndpost = 12000L,   # number of posterior samples
            nskip = 12000L, # number of "warmup" samples to discard
            nchain = 10L,   # number of independent, parallel chains
            ntree = 50L,    # number of trees per chain
            keeptrees = F,
            verbose = F,
            usequants = F,
            k = k_min, base = base_min, power = 2, # hyperparameters tree
            splitprobs = probsmax, keepsampler = F, seed = 22333) 
Samples = fit$yhat.train
dim(Samples)
ids <- sample(1:101, 30, replace = F)
res = list()
for (i in 1:length(ids)) {
  res[[i]]= MCMC_convergence(Samples[,ids[i]],Nchain = 10)
  
}


ypred <- pnorm(fit$yhat.test)
ypred <- colMeans(ypred)

AUC_coBART <- roc(Ytest, ypred, ci = T)
Brier_coBART <- mean((ypred - Ytest)^2)
AUC_coBART
Brier_coBART

##### default BART ####
model <- "flexible"
model <- "rigid"
if (model=="rigid"){k = 1; base = .1; power = 4.0}
if (model=="flexible"){k=2; base = .95; power = 2.0}

fit1 <- bart(x.train = Xtrain, y.train = Ytrain,
            x.test = Xtest,
            ndpost = 12000L,   # number of posterior samples
            nskip = 12000L, # number of "warmup" samples to discard
            nchain = 10L,   # number of independent, parallel chains
            ntree = 50L,    # number of trees per chain
            keeptrees = F,
            verbose = F,
            usequants = F,
            k = k, base = base, power = power, # hyperparameters tree
            splitprobs = c(), keepsampler = F, seed = 3457) 
ypred1 <- pnorm(fit1$yhat.test)
ypred1 <- colMeans(ypred1)
AUC_defBART <- roc(Ytest, ypred1, ci = T)
Brier_defBART <- mean((ypred1 - Ytest)^2)
AUC_defBART
Brier_def = ((ypred1 - Ytest)^2)
Brier_Ebco = ((ypred - Ytest)^2)


wilcox.test(Brier_def,Brier_Ebco, paired = TRUE, alternative = "two.sided")
roc.test(response = Ytest, predictor1=ypred,predictor2=ypred1,
         paired =T)

gc()
#### BART only IPI ####
#model <- "flexible"
model <- "rigid"
if (model=="rigid"){k = 1; base = .1; power = 4.0}
if (model=="flexible"){k=2; base = .95; power = 2.0}

library(dplyr)
Xtrain1 <- select(Xtrain,IPI) 
Xtest1 <- select(Xtest, IPI) 
fit2 <- bart(x.train = Xtrain1, y.train = Ytrain,
             x.test = Xtest1,
             ndpost = 12000L,   # number of posterior samples
             nskip = 12000L, # number of "warmup" samples to discard
             nchain = 10L,   # number of independent, parallel chains
             ntree = 50L,    # number of trees per chain
             keeptrees = F,
             verbose = F,
             usequants = F,
             k = k, base = base, power = power, # hyperparameters tree
             splitprobs = c(), keepsampler = F, seed = 24331) 
ypred2 <- pnorm(fit2$yhat.test)
ypred2 <- colMeans(ypred2)
AUC_defBART <- roc(Ytest, ypred2, ci = T)
Brier_defBART <- mean((ypred2 - Ytest)^2)
AUC_defBART



wilcox.test(Brier_def,Brier_Ebco, paired = TRUE, alternative = "two.sided")
roc.test(response = Ytest, predictor1=ypred,predictor2=ypred2,
         paired =T)

#### ecpc ####
library(ecpc)
load("DatForPaper.Rdata")
Xtrain = dat$Xtrain # feature matrix, 140 features of which 67 copy number variation, 69 mutation, 3 translocation, and the IPI
Ytrain = dat$Ytrain # two year progression free survival (yes = 0, no = 1)
Xtest = dat$Xtest
Ytest = dat$Ytest
CoData = dat$CoData
remove(dat); gc()


CoDatamatrix <- model.matrix(~p.values + Groups,CoData)

fit <- ecpc(Y = Ytrain, X = as.matrix(Xtrain), Z = list(CoDatamatrix),
            maxsel = c(2,5,10,50,80),
            X2 = as.matrix(Xtest), Y2 = Ytest)
a=fit$YpredPost

ecpc_pred <- fit$Ypred
ridge_pred <- unname(fit$Ypredridge[,1])
fit$MSEecpc
fit$MSEridge
fit$beta
fit$intercept

pROC::auc(Ytest,ridge_pred) # ridge
pROC::auc(Ytest,a[,5]) # ecpc

### corf ###
library(CoRF)
library(scam)
library(randomForestSRC)

DFtrain <- data.frame(Ydf=Ytrain,Xdf=Xtrain)
DFtest <- data.frame(Ydf = Ytest, Xdf = Xtest)

## Fitting base rf ##
baseRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",importance=c("none"))

## Fitting co-data model ##
VarsUsed <- baseRF$var.used
#CoDataTrain$pvalsVUmc[is.na(CoDataTrain$pvalsVUmc)] <- mean(CoDataTrain$pvalsVUmc,na.rm=TRUE)

# glm 
CoDataModel  <- glm(VarsUsed/sum(VarsUsed) ~  Groups + p.values, 
                  data=CoData,family='quasibinomial')


# obtaining estimated sampling probabilities
#predswt <- as.numeric(plogis(predict(CoDataModel)))
predswt <- as.numeric(plogis(predict(CoDataModel)))
P <- length(predswt)
predswt2 <- pmax(predswt-1/P,0)
Mtry <- ceiling(sqrt(sum(predswt2!=0)))
predswt2
## Fit new RF with estimated sampling probabilities ##
CoRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",xvar.wt = predswt2, mtry = Mtry)

## Estimating performance for base RF and CoRF ##
preds1 <- predict.rfsrc(baseRF, newdata = DFtest, outcome = "train") # base RF
preds2 <- predict.rfsrc(CoRF, newdata = DFtest, outcome = "train") # CoRF

preds1 <- preds1$predicted # base RF
preds2 <- preds2$predicted # CoRF
#Ytest <- as.numeric(RespValidation)-1
library(pROC)
AUC_baseRF <- roc(Ytest, preds1, ci = T) # base RF
AUC_CoRF <- roc(Ytest, preds2, ci = T) # CoRF
AUC_baseRF
AUC_CoRF
Brier_baseRF <- mean((preds1 - Ytest)^2) # base RF
Brier_CoRF <- (preds2 - Ytest)^2 #CoRF
wilcox.test(Brier_CoRF,Brier_Ebco, paired = TRUE, alternative = "two.sided")
roc.test(response = Ytest, predictor1=ypred,predictor2=preds1,
         paired =T)

##### DART fit #####
library(BART)
set.seed(3456)
fit = pbart(x.train = Xtrain, y.train = Ytrain, x.test = Xtest,
            base = .95, power = 2, k = 2,
            ndpost = 80000L,   # number of posterior samples
            nskip = 20000L, # number of "warmup" samples to discard
            #nchain = 10L,
            keepevery = 20L,
            binaryOffset = mean(Ytrain),
            sparse = T, rho = ncol(Xtrain)/2)

pred_sBART <- fit$prob.test.mean
AUC_1 <- roc(Ytest, pred_sBART, ci = T) # base RF
AUC_1

### Cross-validation for BART ###
#################################
load("DatForPaper.Rdata")
Xtrain = dat$Xtrain # feature matrix, 140 features of which 67 copy number variation, 69 mutation, 3 translocation, and the IPI
Ytrain = dat$Ytrain # two year progression free survival (yes = 0, no = 1)
Xtest = dat$Xtest
Ytest = dat$Ytest

library(dbarts)

cv_trees <- c(50L,150L)
cv_k <- c(1,2,3)
cv_base <- c(0.1,0.5,0.95)


xval <- xbart(Xtrain, Ytrain, n.samples = 5000L, n.test = 5, n.burn = c(1000L, 50L, 50L),
              n.trees = cv_trees,
              k = cv_k,
              power = 2,
              base = c(0.1,0.5,0.95), n.threads = 1L,
              loss = "log",
              method = "k-fold")
save(xval, file = "BartCV_HyperparameterEstimates.Rdata")
load("Results/BartCV_HyperparameterEstimates.Rdata")
xval
xval <- colMeans(xval)
a=which(xval==min(xval), arr.ind = T)
ntree_opt = cv_trees[a[1]]
k_opt = cv_k[a[2]]
base_opt = cv_base[a[3]]
library(dbarts)
fit1 <- bart(x.train = Xtrain, y.train = Ytrain,
             x.test = Xtest,
             ndpost = 100000L,   # number of posterior samples
             nskip = 2000L, # number of "warmup" samples to discard
             nchain = 10L,   # number of independent, parallel chains
             ntree = ntree_opt,    # number of trees per chain
             keeptrees = F,
             verbose = F,
             usequants = F,
             keepevery = 10L, # mcmc thinning,
             binaryOffset = mean(Ytrain),
             k = k_opt, base = base_opt, power = 2, # hyperparameters tree
             splitprobs = c(), keepsampler = F, seed = 5477) 


library(pROC)
ypred_cv <- pnorm(fit1$yhat.test)
ypred_cv <- colMeans(ypred_cv)
AUC_defBART <- roc(Ytest, ypred_cv, ci = T)
Brier_defBART <- mean((ypred_cv - Ytest)^2)
AUC_defBART
load("Results/FinalModelEBCoBARTfit_.95_2_2.Rdata")

library(pROC)
fit = results$FinalFit
Ytest = dat$Ytest
remove(dat)

ypred <- pnorm(fit$yhat.test)
ypred <- colMeans(ypred)
AUC_coBART <- roc(Ytest, ypred, ci = T)
Brier_EBBART <- (ypred - Ytest)^2 # base RF
Brier_CVBART <- (ypred_cv - Ytest)^2 #CoRF

wilcox.test(Brier_EBBART,Brier_CVBART, paired = TRUE, alternative = "two.sided")
roc.test(response = Ytest, predictor1=ypred,predictor2=ypred_cv,
         paired =T)










########## partial dependence plots #########
#############################################
library(dbarts)
library(loo)

load("DatForPaper.Rdata")
Xtrain = dat$Xtrain # feature matrix, 140 features of which 67 copy number variation, 69 mutation, 3 translocation, and the IPI
Ytrain = dat$Ytrain # two year progression free survival (yes = 0, no = 1)
CoData = dat$CoData # Codata matrix, consists of p values (on logit scale) of each feature (Estimated on independent cohort with same outcome), and feature type
Xtest = dat$Xtest
Ytest = dat$Ytest
remove(CoData,dat)

load("Results/FinalModelEBCoBARTResults1_.95_2_2.Rdata")
probsmax <- results$EstProbs
hyps <- results$hyperparameters
k_min = hyps[1]; base_min = hyps[2]; power =hyps[3]
probsmax
pdb_flex <- pdbart(
  x=Xtrain, y = Ytrain,
  xind = c(140),
  ndpost = 50000L,   # number of posterior samples
  nskip = 50000L, # number of "warmup" samples to discard
  nchain = 10L,   # number of independent, parallel chains
  ntree = 50L,    # number of trees per chain
  keeptrees = F,
  verbose = T,
  usequants = F,
  keepevery = 50L, # mcmc thinning
  k = k_min, base = base_min, power = power, # hyperparameters tree
  splitprobs = probsmax, keepsampler = T,
  pl = F, 
)
a_flex=pdb_flex$fd[[1]]

PD_mean = apply(a_flex, MARGIN = 2, mean)
PD_sd = apply(a_flex,MARGIN = 2, sd)

PDS = data.frame(IPI = seq(1,5,1), average = PD_mean, sd = PD_sd, method = rep("flexible EBcoBART",5))

probsmax <- results$EstProbs
hyps <- results$hyperparameters
k_min = hyps[1]; base_min = hyps[2]; power =hyps[3]



pdb_rig <- pdbart(
  x=Xtrain, y = Ytrain,
  xind = c(140),
  ndpost = 50000L,   # number of posterior samples
  nskip = 50000L, # number of "warmup" samples to discard
  nchain = 10L,   # number of independent, parallel chains
  ntree = 50L,    # number of trees per chain
  keeptrees = F,
  verbose = T,
  usequants = F,
  keepevery = 50L, # mcmc thinning
  k = k_min, base = base_min, power = power, # hyperparameters tree
  splitprobs = probsmax, keepsampler = T,
  pl = F, 
)

a_rig=pdb_rig$fd[[1]]

PD_mean_ = apply(a_rig, MARGIN = 2, mean)
PD_sd = apply(a_rig,MARGIN = 2, sd)

PDS_rig = data.frame(IPI = seq(1,5,1), average = PD_mean, sd = PD_sd, method = rep("rigid EBcoBART",5))





pdb1 <- pdbart(
  x=Xtrain, y = Ytrain,
  xind = c(140),
  ndpost = 50000L,   # number of posterior samples
  nskip = 50000L, # number of "warmup" samples to discard
  nchain = 10L,   # number of independent, parallel chains
  ntree = 50L,    # number of trees per chain
  keeptrees = F,
  verbose = T,
  usequants = F,
  keepevery = 50L, # mcmc thinning
  k = 1, base = .1, power = 4, # hyperparameters tree
  splitprobs = c(), keepsampler = T,
  pl = F, 
)

a1=pdb1$fd[[1]]
PD1_mean = apply(a1,MARGIN = 2, mean)
PD1_sd = apply(a1,MARGIN = 2, sd)

PDS1 = data.frame(IPI = seq(1,5,1), average = PD1_mean, sd = PD1_sd, method = rep("rigid defaultBART",5))

pdb2 <- pdbart(
  x=Xtrain, y = Ytrain,
  xind = c(140),
  ndpost = 50000L,   # number of posterior samples
  nskip = 50000L, # number of "warmup" samples to discard
  nchain = 10L,   # number of independent, parallel chains
  ntree = 50L,    # number of trees per chain
  keeptrees = F,
  verbose = T,
  usequants = F,
  keepevery = 50L, # mcmc thinning
  k = 2, base = .95, power = 2, # hyperparameters tree
  splitprobs = c(), keepsampler = T,
  pl = F, 
)

a2=pdb2$fd[[1]]
PD1_mean = apply(a2,MARGIN = 2, mean)
PD1_sd = apply(a2,MARGIN = 2, sd)

PDS2 = data.frame(IPI = seq(1,5,1), average = PD1_mean, sd = PD1_sd, method = rep("flexible defaultBART",5))

PD_tot = rbind.data.frame(PDS,PDS_rig,PDS1,PDS2)

save(PD_tot, file = "PartialDependence_tot_EBcoBART.Rdata")
