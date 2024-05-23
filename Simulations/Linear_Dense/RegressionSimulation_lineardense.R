setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Simulations/Linear_Dense")
library(dbarts)
library(loo)

FiniteSum <- function(x) {
  sum(x[is.finite(x)])
} 
LikelihoodCont <- function(Ypred, Y,sigma){
  # returns the bernoulli likelihoods for the sampled posterior parameters
  # output is a N x m matrix, with N the number of samples and m the number of mc samples
  
  loglik <- -(0.5*(1/sigma^2))*(sweep(Ypred,2,Y)^2)-.5*log(sigma^2)-.5*log(2*pi)
  return(loglik)
}




N = 200
Ntest = 500
p = 500
sigma = 1
#set.seed(1)
#betas1 <- rexp(p, 1)
#betas1 <- sort(betas1,decreasing = T)
#save(betas1, file = "betas_for_sim.Rdata")
load("betas_for_sim.Rdata")


#set.seed(1)
#betas2 = betas1 + rnorm(p, mean = 0, sd = 0.2*sd(betas1))

#save(betas2, file = "CoDataLinDense.Rdata")
load("CoDataLinDense.Rdata")
Codata = data.frame(betas = betas2)
remove(betas2)
library(mvtnorm)
set.seed(12)
Xtest <- rmvnorm(Ntest, mean = rep(0, p), #sigma = CorrX,
                 method= "chol", pre0.9_9994 = FALSE, checkSymmetry = TRUE)
Ytest <- Xtest %*% betas1 + rnorm(Ntest, mean = 0, sd = sigma)

probs <- c(rep(1/p,p))
k =2; base = .95; power = 2.0 # Tree parameters
nu <- 10 ; quant <- .75 # Error variance hyperparameters

nIter <- 10
nrep <-500
GroupProbabilities <- vector("list", length = nrep)
PMSEstest = matrix(NA, nrow = nrep, ncol = nIter)
WAICs = matrix(NA, nrow = nrep, ncol = nIter)


gc()
for (j in 1:nrep){
  print(paste("Sim","Iter",j,sep = " "))
  
  set.seed(j^3+239)
  X <- rmvnorm(N, mean = rep(0, p), #sigma = CorrX,
               method= "chol", pre0.9_9994 = FALSE, checkSymmetry = TRUE)
  Y <- X %*% betas1 + rnorm(N, 0, sigma)
  
  
  # hyperparameter initialization
  probs <- c(rep(1/p,p))
  sigest <- sd(Y)*0.667
  
  # storage containers
  GroupProbs <- matrix(NA, nrow = nIter+1, ncol = p)
  GroupProbs[1,] <- c(rep(1/p,p)) 
  row.names(GroupProbs) <- c(paste("Iter",0:(nIter), sep = " "))
  PMSEtestVector <- c()
  WAICVector <- c()
  
  for (i in 1:nIter) {
    print(paste("iteration",i,sep = " "))
    #print(paste("Probs is:", probs[1], sep = " "))
    fit <- bart(x.train = X, y.train = Y,
                x.test = Xtest,
                ndpost = 12000L,   # number of posterior samples
                nskip = 12000L, # number of "warmup" samples to discard
                nchain = 10L,   # number of independent, parallel chains
                ntree = 50L,    # number of trees per chain
                keeptrees = F,
                verbose = F,
                k = k, base = base, power = power, # hyperparameters tree
                sigest = sigest,sigdf = nu, sigquant = quant,  # hyperparameters error variance
                #proposalprobs = c(0.5,0.5,0,0.5), # relative prob of GROW, PRUNE, CHANGE, and SWAP MH step
                splitprobs = probs) 
    
    # Update Groupprobs
    VarsUsed <- colSums(fit$varcount)
    VarsUsed <- VarsUsed/sum(VarsUsed)
    coDataModel <- glm(VarsUsed ~ 1 + betas,
                       data=Codata,family=quasibinomial)
    #print(coDataModel)
    probs <- predict(coDataModel, type="response", newdata = Codata)
    probs <- unname(probs)
    probs[is.na(probs)] <- 0
    GroupProbs[i+1,] <- probs
    
    PMSEtestVector[i] <- mean((fit$yhat.test.mean-Ytest)^2)
    
    Ypred = fit$yhat.train
    LogLikMatrix = LikelihoodCont(Ypred = Ypred, Y = Y, sigma = fit$sigma)
    WAICVector[i] <- suppressWarnings(waic(LogLikMatrix)$estimates[3,1])
  }
  GroupProbabilities[[j]] <- GroupProbs
  
  PMSEstest[j,] <- PMSEtestVector
  WAICs[j,] <- WAICVector
  
}

name <- paste("LinearDense1",N, k, base, power,"EBBart_WAIC.Rdata", sep = "_")
results <- list(GroupProbs = GroupProbabilities,TestPerformance = PMSEstest, 
                WAIC = WAICs, CoData = Codata)
save(results, file = name)

RatioPMSEtrain <- c()
RatioPMSEtest <- c()
nrep <- length(minima)
var(Ytest)
RatioPMSEtest <- c()
for (i in 1:nrep){
  RatioPMSEtest[i] = PMSEstest[i,1]/PMSEstest[i,minima[i]]
  
}
RatioPMSEtest

######### DART #########
########################
setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Simulations/Linear_Dense")
N = 100
Ntest = 500
p = 500
sigma = 1
load("betas_for_sim.Rdata")

library(mvtnorm)
set.seed(12)
Xtest <- rmvnorm(Ntest, mean = rep(0, p), #sigma = CorrX,
                 method= "chol", pre0.9_9994 = FALSE, checkSymmetry = TRUE)
Ytest <- Xtest %*% betas1 + rnorm(Ntest, mean = 0, sd = sigma)

probs <- c(rep(1/p,p))
k =2; base = .95; power = 2.0 # Tree parameters
nu <- 10 ; quant <- .75 # Error variance hyperparameters
nrep <-500

nchain = 10
theta = ncol(Xtest)

nrep <- 500
library(BART)
PMSEstest <- c()
for (j in 1:nrep){    # for loop representing simulated data sets
  print(paste("Sim","Iter",j,sep = " "))
  
  # simulate training data, either setting 1 or setting 2
  set.seed(j^3+239)
  X <- rmvnorm(N, mean = rep(0, p), #sigma = CorrX,
               method= "chol", pre0.9_9994 = FALSE, checkSymmetry = TRUE)
  Y <- X %*% betas1 + rnorm(N, 0, sigma)
  
  sigest <- sd(Y)*0.667
  
  pred_sBART = matrix(NA,nrow=nchain,ncol = nrow(Xtest))
  for (i in 1:nchain) {
    
    fit = wbart(x.train = X, y.train = Y, x.test = Xtest,
                base = base, power = power, k = k,
                ntree = 50,
                ndpost = 5000L,   # number of posterior samples
                nskip = 5000L, # number of "warmup" samples to discard
                sparse = T, theta = theta,printevery =10000)
    pred_sBART[i,] <- fit$yhat.test.mean
    remove(fit)
  }
  preds <- colMeans(pred_sBART)
  PMSEstest[j] <- mean((preds-Ytest)^2) 
  remove(X,Y,preds)
}
mean(PMSEstest)
name <- paste("LinDense",N, "DART.Rdata", sep = "_")
save(PMSEstest,file = name)






######## ECPC and CoRF ######
#############################

library(ecpc); library(mvtnorm); library(randomForestSRC)
nrep = 500
load("betas_for_sim.Rdata")
load("CoDataLinDense.Rdata")
Z_all = list(matrix(betas2)) # for ecpc
Codata = data.frame(betas=betas2) # for CoRF
p <- 500
sigma <- 1
N <- 100
Ntest <- 500
set.seed(12)

Xtest <- rmvnorm(Ntest, mean = rep(0, p), #sigma = CorrX,
                 method= "chol", pre0.9_9994 = FALSE, checkSymmetry = TRUE)
Ytest <- Xtest %*% betas1 + rnorm(Ntest, mean = 0, sd = sigma)
DFtest <- data.frame(Ydf = Ytest, Xdf = Xtest) # for RF

ecpc_res <- c()
ridge_res <- c()
corf_res <- c()
rf_res <- c()


for (j in 1:nrep){
  print(paste("Sim","Iter",j,sep = " "))
  
  set.seed(j^3+239)
  X <- rmvnorm(N, mean = rep(0, p), #sigma = CorrX,
               method= "chol", pre0.9_9994 = FALSE, checkSymmetry = TRUE)
  Y <- X %*% betas1 + rnorm(N, 0, sigma)
  DFtrain <- data.frame(Ydf=Y,Xdf=X)
  
  ### fit ecpc ###
  res <- ecpc(Y=Y, X = X, Z = Z_all, Y2 = Ytest, X2 = Xtest, postselection = F)
  ecpc_res[j] = res$MSEecpc 
  ridge_res[j] = res$MSEridge
  
  ### fit CoRF ###
  baseRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",
                 importance=c("none"), nodesize = 5, setseed=1)
  ## Fitting co-data model ##
  VarsUsed <- baseRF$var.used
  #CoDataTrain$pvalsVUmc[is.na(CoDataTrain$pvalsVUmc)] <- mean(CoDataTrain$pvalsVUmc,na.rm=TRUE)
  
  # glm 
  CoDataModel  <- glm(VarsUsed/sum(VarsUsed) ~  betas, 
                      data=Codata,family='quasibinomial')
  
  # obtaining estimated sampling probabilities
  #predswt <- as.numeric(plogis(predict(CoDataModel)))
  predswt <- as.numeric(plogis(predict(CoDataModel)))
  predswt2 <- pmax(predswt-1/p,0)
  predswt2 <- predswt2/sum(predswt2)
  Mtry <- sum(predswt2!=0)/3
  ## Fit new RF with estimated sampling probabilities ##
  CoRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",
                xvar.wt = predswt2, mtry = Mtry,
                importance = c("none"), nodesize = 5, setseed = 1)
  preds1 <- predict.rfsrc(baseRF, newdata = DFtest, outcome = "train") # base RF
  preds2 <- predict.rfsrc(CoRF, newdata = DFtest, outcome = "train") # CoRF
  preds1 <- preds1$predicted # base RF
  preds2 <- preds2$predicted # CoRF
  rf_res[j] = mean((preds1-DFtest$Ydf)^2) #base RF
  corf_res[j] = mean((preds2-DFtest$Ydf)^2) #CoRF
}

results <- cbind.data.frame(ridge = ridge_res, ecpc = ecpc_res, rf = rf_res, corf = corf_res)
name = "Results_competitors_N=100_LinDense.Rdata"
save(results, file = name)
colMeans(results)

