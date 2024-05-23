setwd("EB_coBART_paper/Application/Cross-validated Performance")

##############################
### loading relevant libraries ###
library(dbarts)
library(loo)
library(multiridge)




load("DatForPaper.Rdata")
X = dat$Xtrain # feature matrix, 140 features of which 67 copy number variation, 69 mutation, 3 translocation, and the IPI
Y = dat$Ytrain # two year progression free survival (yes = 0, no = 1)
CoData = dat$CoData # Codata matrix, consists of p values (on logit scale) of each feature (Estimated on independent cohort with same outcome), and feature type



source('sidefunctions.R')
set.seed(2)
folds <- CVfolds(Y=Ytrain, balance = T, kfold = 10, nrepeat = 3) ## folds for 5 times repeated 10 fold CV
save(folds, file = "folds.Rdata")


nIter = 18; p = ncol(X)
# initialize model in rigid or flexible tree setting
#model <- "rigid"
model <- "flexible"

EBupdates <- T #set to false if not updating alpha and k

# storage containers for results
GroupProbabilities <- vector("list", length = length(folds))
WAICs = matrix(NA, nrow = length(folds), ncol = nIter)


if (EBupdates == T){
  k_estimates = matrix(NA, nrow = length(folds), ncol = nIter)
  alpha_estimates = matrix(NA, nrow = length(folds), ncol = nIter)
}

for (j in 1:length(folds)) {
  ids = folds[[j]]
  Xtrain = makeModelMatrixFromDataFrame(X[-ids,], drop = F)
  Ytrain <- Y[-ids]
  print(paste("Fold",j,sep = " "))
  EstimatedProbs <- matrix(NA, nrow = nIter+1, ncol = p)
  EstimatedProbs[1,] <- c(rep(1/p,p)) 
  row.names(EstimatedProbs) <- c(paste("Iter",0:(nIter), sep = " "))
  
  # storage containers for single fold
  WAICVector <- c()
  if (EBupdates == T){
    kVector <- c()
    alphaVector <- c()}
    
  
  
  # initialization of EBcoBART for single fold
  probs <- c(rep(1/p,p))
  if (model=="rigid"){k = 1; base = .1; power = 4.0}
  if (model=="flexible"){k=2; base = .95; power = 2.0}
  
  for (i in 1:nIter) {
    print(paste("iteration",i,sep = " "))
    fit <- bart(x.train = Xtrain, y.train = Ytrain,
                ndpost = 12000L,                  # number of posterior samples
                nskip = 12000L,                   # number of "warmup" samples to discard
                nchain = 10L,                     # number of independent, parallel chains
                ntree = 50L,                      # number of trees per chain
                keeptrees = EBupdates,
                usequants = F,
                k = k, base = base, power = power, # hyperparameters tree
                splitprobs = probs,                # prob that variable is chosen for split
                combinechains = T,
                seed = 3*i^2+202+3*i) 
    
    # Codata Model
    VarsUsed <- colSums(fit$varcount)
    VarsUsed <- VarsUsed/sum(VarsUsed)
    Ypred = pnorm(fit$yhat.train)
    
    
    
    gc()
    if (EBupdates == T) {
      # Update leaf node parameter
      trees <- extract(fit,"trees", sampleNums=c(500L:1500L,2500:3500))
      k <- EstimateLeafNode(Trees = trees, ntree = 50)[2]
      kVector[i] <- k
      
      # Update alpha
      trees <- trees[c("tree", "sample", "chain", "n","var","value" )]
      trees$depth <- unname(unlist(by(trees, trees[,c("tree", "sample", "chain")], getDepth)))
      base <- optim(base,LikelihoodTreeStructure,beta = power, Trees=trees, method = 'Brent', lower = .00001, upper = .9999999)$par
      
      alphaVector[i] <- base
      remove(trees)
    }
    
    remove(fit)
    
    coDataModel <- glm(VarsUsed ~ .,
                       data=CoData,family=quasibinomial)
    
    probs <- predict(coDataModel, type="response", newdata = CoData)
    probs <- unname(probs)
    probs[is.na(probs)] <- .0000000000000001
    EstimatedProbs[i+1,] <- probs
    remove(coDataModel,VarsUsed)
    gc()
    
    
    # estimate waic
    Ypred[which(Ypred==0)] <- .0000000000000001
    Ypred[which(Ypred==1)] <- .9999999999999999
    #LikelihoodVector[i] <- Likelihood(preds = colMeans(Ypred), response = Ytrain)
    LogLikMatrix = LikelihoodBin(Ypred = Ypred, Y = Ytrain)
    WAICVector[i] <- waic(LogLikMatrix)[[1]][3,1]
    remove(LogLikMatrix,Ypred)
    
    
  }
  
  WAICs[j,] <- WAICVector
  GroupProbabilities[[j]] <- EstimatedProbs
  if (EBupdates == T){
    k_estimates[j,] <- kVector
    alpha_estimates[j,] <- alphaVector
  }
}

results1 <- list(Probs = GroupProbabilities, alpha = alpha_estimates, k = k_estimates, waic = WAICs, folds = folds)
results1 <- list(Probs = GroupProbabilities, waic = WAICs, folds = folds)
name = "CVResults_EBcoBART_lymphoma_2_2_.95_Treeupdate.Rdata"
name <- paste("CVResults_EBcoBART_lymphoma",model,EBupdates,".Rdata", sep = "_")
save(results1, file = name)


# load estimated hyperparameters for all folds
load("CVResults_EBcoBART_lymphoma_2_2_.95_noTreeupdate.Rdata")
load("CVResults_EBcoBART_lymphoma_2_2_.95_Treeupdate.Rdata")

remove(results1)

### With tree updates ###
folds = results1$folds
GroupProbabilities = results1$Probs
alphas = results1$alpha
ks = results1$k
waics = results1$waic
maxTrain = apply(waics,1, function (x) which.min(x)) # find minimum WAICS for each fold
maxTrain

# determine estimated covariate weights at minimum WAIC
MinGroupProbs1 = matrix(NA, nrow = length(maxTrain), ncol = 140)
for (i in 1:length(maxTrain)){
  MinGroupProbs1[i,] <- GroupProbabilities[[i]][maxTrain[i],]
}

# determine estimated alpha at minimum WAIC
alpha_min <- c()
for (i in 1:length(maxTrain)){
  alpha_min[i] <- alphas[i,maxTrain[i]]
}

# determine estimated k (leaf node hyperparameter) at minimum WAIC
k_min <- c()
for (i in 1:length(maxTrain)){
  k_min[i] <- ks[i,maxTrain[i]]
}

### Without tree updates ###
GroupProbabilities = results$Probs
waics = results$waic
maxTrain = apply(waics,1, function (x) which.min(x))

MinGroupProbs = matrix(NA, nrow = length(maxTrain), ncol = 140)
for (i in 1:length(maxTrain)){
  MinGroupProbs[i,] <- GroupProbabilities[[i]][maxTrain[i],]
}

#### Fitting BART models with estimated hyperparmameters and default BART

library(pROC)
#Without tree update
AUCs_coBART <- c()
Briers_coBART <- c()
Liks_coBART <- c()

#With tree update
AUCs_coBART_Tree <- c()
Briers_coBART_Tree <- c()
Liks_coBART_Tree<- c()

# default BART
AUCs_defBART <- c()
Briers_defBART <- c()
Liks_defBART <- c()
remove(alphas,ks,results,results1,GroupProbabilities,waics)
remove(results1)
gc()
source('sidefunctions.R')
library(dbarts)
for (j in 1:length(folds)) {
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = makeModelMatrixFromDataFrame(X[-ids,], drop = F)
  Ytrain <- Y[-ids]
  Xtest = makeModelMatrixFromDataFrame(X[ids,], drop = F)
  Ytest <- Y[ids]
  
  fit <- bart(x.train = Xtrain, y.train = Ytrain,
              x.test = Xtest,
              ndpost = 12000L,   # number of posterior samples
              nskip = 12000L, # number of "warmup" samples to discard
              nchain = 10L,   # number of independent, parallel chains
              ntree = 50L,    # number of trees per chain
              keeptrees = F,
              verbose = F,
              usequants = F,
              k = k_min[j], base = alpha_min[j], power = 2, # hyperparameters tree
              splitprobs = MinGroupProbs1[j,], keepsampler = F) 
  
  ypred <- pnorm(fit$yhat.test)
  ypred <- colMeans(ypred)
  AUCs_coBART_Tree[j] <- roc(Ytest, ypred, ci = F)$auc
  Briers_coBART_Tree[j] <- mean((ypred - Ytest)^2)
  Liks_coBART_Tree[j] <- Likelihood(preds = ypred, response = Ytest)
  
  remove(fit,ypred)
  
  fit1 <- bart(x.train = Xtrain, y.train = Ytrain,
              x.test = Xtest,
              ndpost = 12000L,   # number of posterior samples
              nskip = 12000L, # number of "warmup" samples to discard
              nchain = 10L,   # number of independent, parallel chains
              ntree = 50L,    # number of trees per chain
              keeptrees = F,
              verbose = F,
              usequants = F,
              k = 2, base = .95, power = 2, # hyperparameters tree
              splitprobs = MinGroupProbs[j,], keepsampler = F) 
  
  ypred1 <- pnorm(fit1$yhat.test)
  ypred1 <- colMeans(ypred1)
  AUCs_coBART[j] <- roc(Ytest, ypred1, ci = F)$auc
  Briers_coBART[j] <- mean((ypred1 - Ytest)^2)
  Liks_coBART[j] <- Likelihood(preds = ypred1, response = Ytest)
  
  remove(fit1,ypred1)

  
  fit2 <- bart(x.train = Xtrain, y.train = Ytrain,
              x.test = Xtest,
              ndpost = 120000L,   # number of posterior samples
              nskip = 12000L, # number of "warmup" samples to discard
              nchain = 10L,   # number of independent, parallel chains
              ntree = 50L,    # number of trees per chain
              keeptrees = F,
              verbose = F,
              usequants = F,
              k = 2, base = .95, power = 2, # hyperparameters tree
              splitprobs = c(), keepsampler = F) 
  # R2 chop primary data, R chop codata, PETAL independent test cohort, 2yrs PFS
  ypred2 <- pnorm(fit2$yhat.test)
  ypred2 <- colMeans(ypred2)
  AUCs_defBART[j] <- roc(Ytest, ypred2, ci = F)$auc
  Briers_defBART[j] <- mean((ypred2 - Ytest)^2)
  Liks_defBART[j] <- Likelihood(preds = ypred2, response = Ytest)
  
  remove(fit2,ypred2)
}
defBART = cbind.data.frame(AUC = AUCs_defBART, Brier = Briers_defBART, Likelihood = Liks_defBART)
coBART_tree = cbind.data.frame(AUC = AUCs_coBART_Tree, Brier = Briers_coBART_Tree, Likelihood = Liks_coBART_Tree)
coBART = cbind.data.frame(AUC = AUCs_coBART, Brier = Briers_coBART, Likelihood = Liks_coBART)

results = list(default = defBART, coEBTree = coBART_tree, coEB = coBART)
name = "CVResults_performances_lymphoma_2_2_.95.Rdata"
save(results,file = name)


load("CVResults_performances_lymphoma_2_2_.95.Rdata")
IPIBART = cbind.data.frame(AUC = AUCs_IPI, Brier = Briers_IPI, Likelihood = Liks_IPI)
results[[4]]<-IPIBART
names(results) <- c("default","coEBTree","coEB","IPI")  
results = append(results, IPIBART)
name = "CVResults_performances_lymphoma_2_2_.95.Rdata"
save(results,file = name)



def = results$default
EBTree = results$coEBTree
EB = results$coEB
IPI = results$IPI
folds = results$folds
library(dplyr)
names(X)<-varnames
X1 <- select(X,IPI)
X1 <- unname(X1)
AUCs_IPI <- c()
Briers_IPI <- c()
Liks_IPI <- c()
remove(results,X)
library(dbarts); library(pROC)
for (j in 1:length(folds)) {
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = X1[-ids,]
  Ytrain <- Y[-ids]
  Xtest = X1[ids,]
  Ytest <- Y[ids]
  
  fit <- bart(x.train = Xtrain, y.train = Ytrain,
              x.test = Xtest,
              ndpost = 12000L,   # number of posterior samples
              nskip = 12000L, # number of "warmup" samples to discard
              nchain = 10L,   # number of independent, parallel chains
              ntree = 50L,    # number of trees per chain
              keeptrees = F,
              verbose = F,
              usequants = F,
              k = 2, base = .95, power = 2, # hyperparameters tree
              splitprobs = c(), keepsampler = F) 
  
  ypred <- pnorm(fit$yhat.test)
  ypred <- colMeans(ypred)
  AUCs_IPI[j] <- roc(Ytest, ypred, ci = F)$auc
  Briers_IPI[j] <- mean((ypred - Ytest)^2)
  Liks_IPI[j] <- Likelihood(preds = ypred, response = Ytest)
  remove(fit,ypred)
}
load("CVResults_performances_lymphoma_2_2_.95.Rdata")
IPIBART = cbind.data.frame(AUC = AUCs_IPI, Brier = Briers_IPI, Likelihood = Liks_IPI)
results[[4]]<-IPIBART
names(results) <- c("default","coEBTree","coEB","IPI")  
results = append(results, IPIBART)
name = "CVResults_performances_lymphoma_2_2_.95.Rdata"
save(results,file = name)


### corf ###
library(CoRF)
library(scam)
library(randomForestSRC)
AUC_baseRF <- c()
Brier_baseRF <- c()
AUC_CoRF <- c()
Brier_CoRF <- c()
library(CoRF)
?CoRF
for (j in 1:length(folds)) {
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = X[-ids,]
  Ytrain <- factor(Y[-ids])
  Xtest = X[ids,]
  Ytest <- factor(Y[ids])
  DFtrain <- data.frame(Ydf=Ytrain,Xdf=Xtrain)
  DFtest <- data.frame(Ydf = Ytest, Xdf = Xtest)
  ### Fitting Corf ###
  baseRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",importance=c("none"))
  preds1 <- predict.rfsrc(baseRF, newdata = DFtest, outcome = "train") # base RF
  preds1 <- preds1$predicted[,2] # base RF
  AUC_baseRF[j] <- roc(Ytest, preds1, ci = T)$auc # base RF
  Brier_baseRF[j] <- mean((preds1 - (as.numeric(Ytest)-1))^2) # base RF
  
  
  ## Fitting co-data model ##
  VarsUsed <- baseRF$var.used
  remove(preds1,baseRF)
  CoDataModel <-  glm(VarsUsed/sum(VarsUsed) ~ 1 + Groups + p.values, 
                      data=CoData,family='quasibinomial')
  predswt <- as.numeric(plogis(predict(CoDataModel)))
  P <- length(predswt)
  predswt2 <- pmax(predswt-1/P,0)
  Mtry <- ceiling(sqrt(sum(predswt2!=0)))
  CoRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",xvar.wt = predswt2, mtry = Mtry)
  
  ## Estimating performance for base RF and CoRF ##
  
  preds2 <- predict.rfsrc(CoRF, newdata = DFtest, outcome = "train") # CoRF
  preds2 <- preds2$predicted[,2] # CoRF
  AUC_CoRF[j] <- roc(Ytest, preds2, ci = T)$auc # CoRF
  Brier_CoRF[j] <- mean((preds2 - (as.numeric(Ytest)-1))^2) #CoRF
  remove(preds2,CoRF)
}
baseRF = cbind.data.frame(AUC = AUC_baseRF, Brier = Brier_baseRF)
CoRF = cbind.data.frame(AUC = AUC_CoRF, Brier = Brier_CoRF)

load("CVResults_performances_lymphoma_2_2_.95.Rdata")

results[[5]]<-baseRF
results[[6]]<-CoRF
names(results) <- c("default","coEBTree","coEB","IPI","baseRF","CoRF")  
name = "CVResults_performances_lymphoma_2_2_.95.Rdata"
save(results,file = name)

### ecpc ###
CoDatamatrix <- model.matrix(~p.values + Groups,CoData)
library(ecpc)
?ecpc
AUC_ridge <- c()
AUC_ecpc <- c()
Brier_ridge <- c()
Brier_ecpc <- c()


for (j in 1:length(folds)) {
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = X[-ids,]
  Ytrain <- Y[-ids]
  Xtest = X[ids,]
  Ytest <- Y[ids]
  
  fit <- ecpc(Y = Ytrain, X = as.matrix(Xtrain), Z = list(CoDatamatrix), postselection = F,
              X2 = as.matrix(Xtest), Y2 = Ytest)
  
  print(fit$gamma)
  ypredRidge = fit$Ypredridge[,1]
  AUC_ridge[j] <- roc(Ytest, ypredRidge, ci = T)$auc
  Brier_ridge[j] <- mean((ypredRidge - Ytest)^2)
  remove(ypredRidge)
  
  
  ypredECPC = fit$Ypred
  AUC_ecpc[j] <- roc(Ytest, ypredECPC, ci = T)$auc
  Brier_ecpc[j] <- mean((ypredECPC - Ytest)^2)
  remove(ypredECPC,fit,Xtrain,Ytrain,Xtest,Ytest)
}


Ridge = cbind.data.frame(AUC = AUC_ridge, Brier = Brier_ridge)
ECPC = cbind.data.frame(AUC = AUC_ecpc, Brier = Brier_ecpc)


load("CVResults_performances_lymphoma_2_2_.95.Rdata")

results[[7]]<-Ridge
results[[8]]<-ECPC
names(results) <- c("default","coEBTree","coEB","IPI",
                    "baseRF","CoRF","Ridge","ecpc")  
name = "CVResults_performances_lymphoma_2_2_.95.Rdata"
save(results,file = name)







##### DART fit #####
load("DatForPaper.Rdata")
X = dat$Xtrain # feature matrix, 140 features of which 67 copy number variation, 69 mutation, 3 translocation, and the IPI
Y = dat$Ytrain #

load("CVResults_EBcoBART_lymphoma_2_2_.95_Treeupdate.Rdata")

folds = results$folds; remove(results); gc()




library(BART)

nchain = 10
theta = ncol(X)
AUC <- c()
Brier <- c()
for (j in 1:length(folds)) {
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = X[-ids,]
  Ytrain <- Y[-ids]
  Xtest = X[ids,]
  Ytest <- Y[ids]
  pred_sBART = matrix(NA,nrow=nchain,ncol = length(Ytest))
  for (i in 1:nchain) {
    
    print(paste("Chain",i,sep = " "))
    set.seed(34*i)
    fit = pbart(x.train = Xtrain, y.train = Ytrain, x.test = Xtest,
                base = .95, power = 2, k = 2,
                ndpost = 12000L,   # number of posterior samples
                nskip = 12000L, # number of "warmup" samples to discard
                sparse = T, theta = theta,printevery=50000)
    pred_sBART[i,] <- fit$prob.test.mean
  }
  
  pred_sBART1 = colMeans(pred_sBART)
  
  AUC[j] <- roc(Ytest, pred_sBART1, ci = T)$auc 
  Brier[j] = mean((pred_sBART1-Ytest)^2)
  remove(Xtrain,Ytrain,Xtest,Ytest,ids,pred_sBART,pred_sBART1)
}
mean(AUC)
mean(Brier)
DART = cbind.data.frame(AUC = AUC, Brier = Brier)

load("CVedPerformances_lymphoma_all.Rdata")  
library(rlist)

results = list.append(results2,"DART"= DART)
name = "CVedPerformances_lymphoma_all.Rdata"
save(results,file = name) 
load("CVedPerformances_lymphoma_all.Rdata")  
  