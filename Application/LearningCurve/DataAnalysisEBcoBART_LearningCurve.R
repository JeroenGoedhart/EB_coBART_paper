###### Learning Curves for BART #####

setwd("EB_coBART_paper/Application/LearningCurve")



load("DatForPaper.Rdata")
X = dat$Xtrain # feature matrix, 140 features of which 67 copy number variation, 69 mutation, 3 translocation, and the IPI
Y = dat$Ytrain # two year progression free survival (yes = 0, no = 1)
CoData = dat$CoData # Codata matrix, consists of p values (on logit scale) of each feature (Estimated on independent cohort with same outcome), and feature type
remove(dat)

source('sidefunctions.R') # load side functions consists of
# EstimateLeafNode: empirical bayes estimates of k hyperparameter (leaf node) of BART
# FiniteSum: only sum observations that have finite values
# getDepth: computes the depth of each node of all trees of a BART fit
# LikelihoodBin: computes likelihood values for binary variables based on bernoulli distribution
# LikelihoodTreestructure: compute the prior probability on log scale of a tree structure to occur as a function of
# alpha and beta hyperparameters

Subs <- function(Y,model="logistic",balance=TRUE,ntrain,fixedsubs=TRUE,nrepeat=1){ #response is required for balanced CV
  #response: response vector, length n
  #model: "logistic", "cox", etc
  #balance: should the splits balance levels of the response?
  #kfold: scalar, the number of folds
  #fixedsubs: should the folds be fixed? (for reproducibility)
  #nrepeat: number of repeats of the CV
  #Output: list object with kfold elements containing the sample indices of the left-out samples per fold
  
  #Y <- rbinom(100,1,.95);model="logistic";balance=FALSE;ntrain=90;fixedsubs=TRUE;nrepeat=1
  
  response <- Y
  if(model=="linear") balance <- FALSE
  subrep <- function(rep){
    nsam <- length(response)
    ntest <- nsam-ntrain
    if (fixedsubs) set.seed(3534+rep-1) else set.seed(NULL)
    if (!balance) {
      subs <- sort(sample(1:nsam,ntest))
    }
    else {
      if (model == "logistic"){
        if (class(response) == "factor")
          nev <- which((as.numeric(response) - 1) == 1)
        else nev <- which(response == 1)
      }
      if (model == "cox") nev <- which(response[, 2] == 1)
      nonev <- (1:nsam)[-nev]
      nsamev <- length(nev)
      ntestev <- min(max(1,round(ntest/nsam * nsamev)),(nsamev-1))
      ntestnonev <- ntest - ntestev
      subs <- sort(c(sample(nev,ntestev), sample(nonev,ntestnonev)))
    }
    return(subs)
  }
  return(lapply(1:nrepeat,subrep))
}
nmin = 20; nev = 10
nseq <- round(seq(nmin,nrow(X)-10,length.out=nev))
Subsamples <- vector("list", length = nev)
set.seed(1)
for (i in 1:nev) {
 Subsamples[[i]] <- Subs(Y = Y, model = "logistic", 
                         balance = T, ntrain = nseq[i],
                         nrepeat = 1)[[1]]
  
}

gc()
nIter = 15; p = ncol(X)
# storage containers
GroupProbabilities <- vector("list", length = nev)
WAICs = matrix(NA, nrow = nev, ncol = nIter)
k_estimates = matrix(NA, nrow = nev, ncol = nIter)
alpha_estimates = matrix(NA, nrow = nev, ncol = nIter)

model <- "flexible"
model <- "rigid"

gc()
library(dbarts); library(loo)

for (j in 1:nev) {
  ids = Subsamples[[j]]
  Xtrain = makeModelMatrixFromDataFrame(X[-ids,], drop = F)
  Ytrain <- Y[-ids]
  print(paste("Sample size",length(Y)-length(ids),sep = " "))
  EstimatedProbs <- matrix(NA, nrow = nIter+1, ncol = p)
  EstimatedProbs[1,] <- c(rep(1/p,p)) 
  row.names(EstimatedProbs) <- c(paste("Iter",0:(nIter), sep = " "))
  
  WAICVector <- c()
  kVector <- c()
  alphaVector <- c()
  
  probs <- c(rep(1/p,p))
  if (model=="rigid"){k = 1; base = .1; power = 4.0}
  if (model=="flexible"){k=2; base = .95; power = 2.0}
  
  
  for (i in 1:nIter) {
    print(paste("iteration",i,sep = " "))
    fit <- bart(x.train = Xtrain, y.train = Ytrain,
                ndpost = 12000L,                   # number of posterior samples
                nskip = 12000L,                    # number of "warmup" samples to discard
                nchain = 10L,                      # number of independent, parallel chains
                ntree = 50L,                       # number of trees per chain
                keeptrees = T,
                verbose = F,
                usequants = F,
                k = k, base = base, power = power, # hyperparameters tree
                splitprobs = probs,                # prob that variable is chosen for split
                combinechains = T,
                #proposalprobs = c(0.5,0.25,0.25,0.5),
                seed = 3*i^2+202+3*i) 
    
    # Codata Model
    VarsUsed <- colSums(fit$varcount)
    VarsUsed <- VarsUsed/sum(VarsUsed)
    Ypred = pnorm(fit$yhat.train)
    trees <- extract(fit,"trees",chainNum = c(1:10), sampleNum=c(1:100,1000:1100,1900:2000))
    
    remove(fit)
    gc()
    
    # Update leaf node parameter
    k <- EstimateLeafNode(Trees = trees, ntree = 50)[2]
    kVector[i] <- k
    
    # Update alpha
    trees <- trees[c("tree", "sample", "chain", "n","var","value" )]
    trees$depth <- unname(unlist(by(trees, trees[,c("tree", "sample", "chain")], getDepth)))
    base <- optim(base,LikelihoodTreeStructure,beta = power, Trees=trees, method = 'Brent', lower = .00001, upper = .9999999)$par
    
    alphaVector[i] <- base
    remove(trees)
    
    
    
    
    coDataModel <- glm(VarsUsed ~ .,
                       data=CoData,family=quasibinomial)
    #print(coDataModel)
    probs <- predict(coDataModel, type="response", newdata = CoData)
    probs <- unname(probs)
    probs[is.na(probs)] <- .0000000000000001
    EstimatedProbs[i+1,] <- probs
    remove(coDataModel,VarsUsed)
    gc()
    # estimate predictive performance
    #Ypred <- colMeans(pnorm(fit$yhat.test))
    #AUCtestVector[i] <- roc(Ytest,colMeans(pnorm(fit$yhat.test)), quiet = T)$auc
    #BrierScoreVector[i] <- mean((Ypred-Ytest)^2)
    
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
  k_estimates[j,] <- kVector
  alpha_estimates[j,] <- alphaVector
}

results<- list(Probs = GroupProbabilities, alpha = alpha_estimates, k = k_estimates, waic = WAICs, Subsamples = Subsamples)
name = paste("LearningCurveResults_EBcoBART_EstHyperpars", model,".Rdata",sep = "_")
save(results, file = name)
load("DatForPaper.Rdata")
X = dat$Xtrain
Y = dat$Ytrain

remove(dat); gc()
load("LearningCurveResults_EBcoBART_EstHyperpars_1_4_.1_Treeupdate.Rdata")
GroupProbabilities = results$Probs
alphas = results$alpha
ks = results$k
waics = results$waic

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

maxTrain = apply(waics,1, function (x) minimum(x))
maxTrain

MinGroupProbs = matrix(NA, nrow = length(maxTrain), ncol = 140)
for (i in 1:length(maxTrain)){
  MinGroupProbs[i,] <- GroupProbabilities[[i]][maxTrain[i],]
}

alpha_min <- c()
for (i in 1:length(maxTrain)){
  alpha_min[i] <- alphas[i,maxTrain[i]]
}
k_min <- c()
for (i in 1:length(maxTrain)){
  k_min[i] <- ks[i,maxTrain[i]]
}
alpha_min
k_min
remove(maxTrain)
remove(results,GroupProbabilities, alphas, ks, waics,i)

nmin = 20; nev = 10; nrep = 50
nseq <- round(seq(nmin,91,length.out=nev))
row.names(MinGroupProbs) <- nseq 
names(alpha_min) <- nseq
names(k_min) <- nseq


p = nrow(X)
# storage
AUCs_EBCo <- vector("list", length = nev)
Briers_EBCo <- vector("list", length = nev) 
AUCs_def <- vector("list", length = nev) 
Briers_def <- vector("list", length = nev)
model <- "rigid"
#model <- "flexible"

library(dbarts); library(pROC)

for (i in 1:nev) {
  print(paste("Sample size",nseq[i],sep = " "))
  set.seed(i^3+5)
  Subsamps <- Subs(Y = Y, model = "logistic", 
                   balance = T, ntrain = nseq[i],
                   nrepeat = nrep)
  alpha = alpha_min[i]; k = k_min[i]; 
  if (model=="rigid"){power = 4.0}
  if (model=="flexible"){power = 2.0}
  
  probs = MinGroupProbs[i,]
  AUC_nev_EBCo <- c()
  Brier_nev_EBCo <- c()
  AUC_nev_def <- c()
  Brier_nev_def <- c()
  for (j in 1:nrep){
    ids <- Subsamps[[j]]
    Xtrain = makeModelMatrixFromDataFrame(X[-ids,], drop = F)
    Xtest = makeModelMatrixFromDataFrame(X[ids,], drop = F)
    Ytrain <- Y[-ids]; Ytest <- Y[ids]
    print(paste("repeat",j,sep = " "))
    #EBcoBART
    fit <- bart(x.train = Xtrain, y.train = Ytrain,
                x.test = Xtest,
                ndpost = 10000L,   # number of posterior samples
                nskip = 2000L, # number of "warmup" samples to discard
                nchain = 10L,   # number of independent, parallel chains
                ntree = 50L,    # number of trees per chain
                keeptrees = F,
                verbose = F,
                usequants = F,
                #keepevery = 10L, # mcmc thinning,
                #binaryOffset = mean(Ytrain),
                k = k, base = alpha, power = power, # hyperparameters tree
                splitprobs = probs, keepsampler = F) 
    
    ypred <- pnorm(fit$yhat.test)
    ypred <- colMeans(ypred)
    AUC_nev_EBCo[j] <- roc(Ytest, ypred, ci = F)$auc
    Brier_nev_EBCo[j] <- mean((ypred - Ytest)^2)
    remove(fit,ypred)
    # default BART
    fit1 <- bart(x.train = Xtrain, y.train = Ytrain,
                x.test = Xtest,
                ndpost = 10000L,   # number of posterior samples
                nskip = 2000L, # number of "warmup" samples to discard
                nchain = 10L,   # number of independent, parallel chains
                ntree = 50L,    # number of trees per chain
                keeptrees = F,
                verbose = F,
                usequants = F,
                k = 2, base = .95, power = power, # hyperparameters tree
                splitprobs = c(), keepsampler = F) 
    
    ypred1 <- pnorm(fit1$yhat.test)
    ypred1 <- colMeans(ypred1)
    AUC_nev_def[j] <- roc(Ytest, ypred1, ci = F)$auc
    Brier_nev_def[j] <- mean((ypred1 - Ytest)^2)
    print(c(AUC_nev_EBCo[j],AUC_nev_def[j]))
    remove(fit1,ypred1,Xtrain,Xtest,Ytrain,Ytest)
    gc()
    
  }
  #EBcoBART
  AUCs_EBCo[[i]] <- AUC_nev_EBCo
  Briers_EBCo[[i]] <- Brier_nev_EBCo
  # defBART
  AUCs_def[[i]] <- AUC_nev_def 
  Briers_def[[i]] <- Brier_nev_def
  
  remove(AUC_nev_EBCo,Brier_nev_def,AUC_nev_def,Brier_nev_EBCo); gc()
  
}
names(AUCs_EBCo) <- nseq; names(Briers_EBCo) <- nseq
names(AUCs_def) <- nseq; names(Briers_def) <- nseq
AUCs_EBCo
AUC_nev_EBCo[1]
res = cbind(alpha_min,k_min)

results1 <- list(AUCs_EBCo = AUCs_EBCo, AUCs_def = AUCs_def, 
                 Briers_EBCo = Briers_EBCo, Briers_def = Briers_def)
name = paste("LearningCurveResults_EBcoBART_Performances",model,".Rdata",sep = "_")
save(results1, file = name)
getwd()
load("LearningCurveResults_EBcoBART_Performances_1_4_.1_Treeupdate.Rdata")
Briers_def <- results1$Briers_def
Briers_EBco <- results1$Briers_EBCo
AUCs_def <- results1$AUCs_def
AUCs_EBco <- results1$AUCs_EBCo
Brier_def <- unlist(lapply(Briers_def, mean))
Brier_EBco <- unlist(lapply(Briers_EBco, mean))
AUC_def <- unlist(lapply(AUCs_def, mean))
AUC_EBco <- unlist(lapply(AUCs_EBCo, mean))



