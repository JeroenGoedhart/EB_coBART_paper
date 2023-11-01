#### comparison with cross-validation ####

setwd("C:/Users/VNOB-0732/Desktop/R files/Bayesian Additive Regression Trees/Grouping Information/GroupingBART/SimulationsEBBart")

FiniteSum <- function(x) {
  sum(x[is.finite(x)])
} 
LikelihoodCont <- function(Ypred, Y,sigma){
  # returns the bernoulli likelihoods for the sampled posterior parameters
  # output is a N x m matrix, with N the number of samples and m the number of mc samples
  
  loglik <- -(0.5*(1/sigma^2))*(sweep(Ypred,2,Y)^2)-.5*log(sigma^2)-.5*log(2*pi)
  return(loglik)
}
GroupProb <- function(VarCounts,Str_Group){
  N_G <- length(Str_Group)
  Probs_G <- c()
  Probs_G[1] <-sum(VarCounts[1:Str_Group[1]])/sum(VarCounts)
  for (i in 2:N_G) {
    sum_min <- sum(Str_Group[1:(i-1)])+1
    sum_max <- sum(Str_Group[1:i])
    #print(c(sum_min,sum_max))
    #print(c(sum_min,sum_max))
    Probs_G[i]<- sum(VarCounts[sum_min:sum_max])/sum(VarCounts)
  }
  ProbUpdate <- c()
  for (i in 1:N_G) {
    ProbUpdate <- append(ProbUpdate,values = rep(Probs_G[i]/Str_Group[i], Str_Group[i]))
    
  }
  return(list(Probs_G/sum(Probs_G),ProbUpdate))
}



library(dbarts)
library(loo)
library(multiridge)
library(mvtnorm)

## Set-up ##
p <- 500
sigma <- 1.0
N <- 100
G <- 5   #number of groups
GroupStructure <- c(rep(p/G,G))
Ntest =500


### Data setting 1: sparse and nonlinear ###
g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,101] - 0.5)^2 +
    10 * x[,102] + 10 * x[,3] #+ 10 * sin(pi * x[,201] * x[,202]) + 20 * (x[,301] - 0.5)^2 +
  #+ 10 * x[,203] + 10 * x[,302]
}



Xtest1 <- matrix(runif(Ntest*p),Ntest,p)
Ytest1 <- g(Xtest1) + rnorm(Ntest,0,sigma)


# Setting up simulation
probs <- c(rep(1/p,p))
k =2; base = .95; power = 2.0 # Tree parameters
nu <- 10 ; quant <- .75 # Error variance hyperparameters
nchain = 10

nIter <- 15
nrep <-100
GroupProbabilities <- vector("list", length = nrep)
#Treestats <- vector("list", length = nrep)
PMSEsCV = matrix(NA, nrow = nrep, ncol = nIter)
PMSEstest = matrix(NA, nrow = nrep, ncol = nIter)
#PMSEstrain = matrix(NA, nrow = nrep, ncol = nIter)
#kupdates = matrix(NA, nrow = nrep, ncol = nIter+1)
WAICs = matrix(NA, nrow = nrep, ncol = nIter)


gc()
for (j in 1:nrep){
  print(paste("Sim","Iter",j,sep = " "))
  
  set.seed(j^3+239)
  #X <- rmvnorm(N, mean = rep(0, p), #sigma = CorrX,
   #            method= "chol", pre0.9_9994 = FALSE, checkSymmetry = TRUE)
  X <- matrix(runif(N * p), N, p)
  #Y <- X %*% betas1 + rnorm(N, 0, sigma)
  Y <- g(X)+ rnorm(N, 0, sigma)
  folds <- CVfolds(Y, kfold = 5, nrepeat = 1)
  # hyperparameter initialization
  probs <- c(rep(1/p,p))
  #k <- 2
  #base = .95 ; power = 2
  #print(c(probs[1],k,base,power))
  sigest <- sd(Y)*0.667
  
  # storage containers
  GroupProbs <- matrix(NA, nrow = nIter+1, ncol = G)
  GroupProbs[1,] <- c(rep(1/G,G)) 
  row.names(GroupProbs) <- c(paste("Iter",0:(nIter), sep = " "))
  PMSEtestVector <- c()
  PMSECV_Vector <- c()
  WAICVector <- c()
  #k_Update <- c(k)
  #Treepars <- c(base,power)
  
  
  
  #### no cv loop ####
  for (i in 1:nIter) {
    print(paste("iteration",i,sep = " "))
    #print(paste("Probs is:", probs[1], sep = " "))
    fit <- bart(x.train = X, y.train = Y,
                x.test = Xtest1,
                ndpost = 10000L,   # number of posterior samples
                nskip = 2000L, # number of "warmup" samples to discard
                nchain = nchain,   # number of independent, parallel chains
                ntree = 50L,    # number of trees per chain
                keeptrees = F,
                verbose = F,
                keepevery = 2L, # mcmc thinning
                k = k, base = base, power = power, # hyperparameters tree
                sigest = sigest,sigdf = nu, sigquant = quant,  # hyperparameters error variance
                #proposalprobs = c(0.5,0.5,0,0.5), # relative prob of GROW, PRUNE, CHANGE, and SWAP MH step
                splitprobs = probs) 
    
    # Update Groupprobs
    a<-GroupProb(VarCounts = colSums(fit$varcount), Str_Group = GroupStructure)
    probs <- a[[2]]
    GroupProbs[i+1,] <- a[[1]]
    
    # estimate WAIC
    Ypred = fit$yhat.train
    LogLikMatrix = LikelihoodCont(Ypred = Ypred, Y = Y, sigma = fit$sigma)
    WAICVector[i] <- suppressWarnings(waic(LogLikMatrix)$estimates[3,1])
    
    # estimate test performance
    PMSEtestVector[i] =mean((fit$yhat.test.mean-Ytest1)^2)
    
    # Estimate cv performance
    PMSEfolds <- c()
    print("CV estimation")
    for (l in 1:length(folds)) {
      print(paste("fold",l,sep = " "))
      ids <- folds[[l]]
      Xtrain = X[-ids,]; Ytrain = Y[-ids]
      Xtest = X[ids,]; Ytest = Y[ids]
      fit1 <- bart(x.train = Xtrain, y.train = Ytrain,
                  x.test = Xtest,
                  ndpost = 10000L,   # number of posterior samples
                  nskip = 2000L, # number of "warmup" samples to discard
                  nchain = nchain,   # number of independent, parallel chains
                  ntree = 50L,    # number of trees per chain
                  keeptrees = F,
                  verbose = F,
                  keepevery = 2L, # mcmc thinning
                  k = k, base = base, power = power, # hyperparameters tree
                  sigest = sigest,sigdf = nu, sigquant = quant,  # hyperparameters error variance
                  #proposalprobs = c(0.5,0.5,0,0.5), # relative prob of GROW, PRUNE, CHANGE, and SWAP MH step
                  splitprobs = probs)
      PMSEfolds[l] =mean((fit1$yhat.test.mean-Ytest)^2)
      remove(fit1)
      remove(Xtrain,Ytrain,Xtest,Ytest)
    }
    PMSECV_Vector[i] <- mean(PMSEfolds)
  }

  GroupProbabilities[[j]] <- GroupProbs
  
  PMSEstest[j,]<- PMSEtestVector
  PMSEsCV[j,] <- PMSECV_Vector
  WAICs[j,] <- WAICVector
  
  remove(GroupProbs,probs,PMSEtestVector,WAICVector,PMSECV_Vector)
  gc()
}
name <- paste("Friedman",G,N, k, base, power,"CVversusWAIC.Rdata", sep = "_")
results <- list(GroupProbs = GroupProbabilities,TestPerformance = PMSEstest, 
                WAIC = WAICs,CVPerf = PMSEsCV)
save(results, file = name)

CV <- results$CVPerf
Test <- results$TestPerformance
WAIC <- results$WAIC
GroupProbs <- results$GroupProbs
remove(results)

minCV <-apply(CV,1, function (x) which(x == min(x)))

minWAIC <-apply(WAIC,1, function (x) which(x == min(x)))

minTest <- apply(Test,1, function (x) which(x == min(x)))

RatioPMSEtest <- c()

nrep <- length(minCV)
for (i in 1:nrep){
  RatioPMSEtest[i] = Test[i,minWAIC[i]]/Test[i,minCV[i]]
}


col = "#395D9CFF"
col= "#0B0405FF"
col= "#60CEACFF"

RatioPMSEtest
which(RatioPMSEtest>=1.03)


name <- paste("RatioPerf_WAICvsCV.pdf", sep = "_")
pdf(name, width=7,height=4)

boxplot(RatioPMSEtest,
        col = "#00A087B2", 
        outline = F, font.lab = 2, xaxt="n", yaxt="n")
axis(2, family = "serif", font.lab=2, font=2)
#axis(1, cex = 0.1, family = "serif", font.lab=2, font=2)
title(ylab = "PMSE Ratio", xlab = "", family = "serif", font.lab=2, font=2, cex.lab=1.1)
dev.off()
nrep = 100; ncol = 5

MinGroupProbsCV = matrix(NA, nrow = nrep, ncol = ncol)
for (i in 1:nrep){
  MinGroupProbsCV[i,] <- GroupProbs[[i]][minCV[i],]
}

GroupProbs[[1]]


MinGroupProbsWAIC = matrix(NA, nrow = nrep, ncol = ncol)
for (i in 1:nrep){
  MinGroupProbsWAIC[i,] <- GroupProbs[[i]][minWAIC[i],]
}

diff = MinGroupProbsWAIC-MinGroupProbsCV

max(diff)
min(diff)
length(which(diff[,1]==0))
library(viridis)
cols = viridis(5, direction = -1, option = "H")
names<- c(paste("G",1:5, sep = ""))
colMeans(diff)
name <- paste("GroupProbDiff_WAICvsCV.pdf", sep = "_")
pdf(name, width=7,height=4)
boxplot(diff, col = cols, ylim = c(0,1), yaxt = "n", xaxt = "n", outline = F)
axis(2, family = "serif", font.lab=2, font=2)
axis(1, family = "serif", font.lab=2, font=2, labels = names, at = seq(1,5,1))
title(ylab = "Group Prob. Diff.", xlab = "", family = "serif", font.lab=2, font=2, cex.lab=1.1)

dev.off()



boxplot(RatioPMSEtest,outline = F)
mean(RatioPMSEtest)


mean(minWAIC-minCV)

warnings()
