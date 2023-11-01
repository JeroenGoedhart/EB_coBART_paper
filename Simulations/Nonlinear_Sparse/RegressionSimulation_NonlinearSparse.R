gc()

setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Simulations/Nonlinear_Sparse")

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


## Set-up ##
p <- 500
sigma <- 1.0
N <- 100
G <- 5   #number of groups
GroupStructure <- c(rep(p/G,G))
Ntest =1000


### Data setting 1: sparse and nonlinear ###
g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,101] - 0.5)^2 + 10 * x[,102] + 10 * x[,3]
}

g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 10 * x[,5] +
  10 * sin(pi * x[,101] * x[,102]) + 20 * (x[,103] - 0.5)^2 + 10 * x[,104] + 10 * x[,105] +
  10 * sin(pi * x[,201] * x[,202]) + 20 * (x[,203] - 0.5)^2 + 10 * x[,204] + 10 * x[,205] +
  10 * sin(pi * x[,301] * x[,302]) + 20 * (x[,303] - 0.5)^2 + 10 * x[,304] + 10 * x[,305] +
  10 * sin(pi * x[,401] * x[,402]) + 20 * (x[,403] - 0.5)^2 + 10 * x[,404] + 10 * x[,405]
}


Xtest <- matrix(runif(Ntest*p),Ntest,p)
Ytest <- g(Xtest) + rnorm(Ntest,0,sigma)


source('SideFunctions_EmpiricalBayes_BART.R')



# Setting up simulation
probs <- c(rep(1/p,p))
k =2; base = .95; power = 2.0 # Tree parameters
nu <- 10 ; quant <- .75 # Error variance hyperparameters

nIter <- 10
nrep <- 500

## Result containers 
GroupProbabilities <- vector("list", length = nrep)
PMSEstest = matrix(NA, nrow = nrep, ncol = nIter)
WAICs = matrix(NA, nrow = nrep, ncol = nIter)

for (j in 1:nrep){    # for loop representing simulated data sets
  print(paste("Sim","Iter",j,sep = " "))
  
  # simulate training data, either setting 1 or setting 2
  set.seed(j^3+239)
  X <- matrix(runif(N * p), N, p)
  Y <- g(X)+ rnorm(N, 0, sigma)
  
  # hyperparameter initialization
  probs <- c(rep(1/p,p))
  sigest <- sd(Y)*0.667
  
  # storage containers
  GroupProbs <- matrix(NA, nrow = nIter+1, ncol = G)
  GroupProbs[1,] <- c(rep(1/G,G)) 
  row.names(GroupProbs) <- c(paste("Iter",0:(nIter), sep = " "))
  PMSEtestVector <- c()
  WAICVector <- c()
  
  
  #### Entering EBcoBART (iteratively updating splitprobs hyperparameter of dbarts) ####
  for (i in 1:nIter) {
    print(paste("iteration",i,sep = " "))
    fit <- bart(x.train = X, y.train = Y, # training data
                x.test = Xtest,  # testing data
                ndpost = 120000L,   # number of posterior samples
                nskip = 12000L, # number of "warmup" samples to discard
                nchain = 10L,   # number of independent, parallel chains
                ntree = 50L,    # number of trees per chain
                keeptrees = F,
                verbose = F,
                k = k, base = base, power = power, # hyperparameters tree
                sigest = sigest,sigdf = nu, sigquant = quant,  # hyperparameters error variance
                splitprobs = probs   # hyperparameter that will be updated using EB and co-data
                ) 
    
    # Update Groupprobs
    a<-GroupProb(VarCounts = colSums(fit$varcount), Str_Group = GroupStructure) # estimate the group probability
    probs <- a[[2]]
    GroupProbs[i+1,] <- a[[1]]
    print(a[[1]][1])
    
    # Estimate  test performance
    PMSEtestVector[i] <- mean((fit$yhat.test.mean-Ytest)^2)
    
    # Estimate WAIC
    Ypred = fit$yhat.train
    LogLikMatrix = LikelihoodCont(Ypred = Ypred, Y = Y, sigma = fit$sigma)
    WAICVector[i] <- suppressWarnings(waic(LogLikMatrix)$estimates[3,1])
    
  }
  GroupProbabilities[[j]] <- GroupProbs
  
  
  PMSEstest[j,] <- PMSEtestVector
  WAICs[j,] <- WAICVector
  remove(GroupProbs,probs,WAICVector)
  gc()
}
 
name <- paste("Friedman",G,N, k, base, power,"EBBart_WAIC.Rdata", sep = "_")
results <- list(GroupProbs = GroupProbabilities,TestPerformance = PMSEstest, waic = WAICs)

save(results, file = name)


#############################################################
########################### ecpc ############################
#############################################################
## Set-up ##
p <- 500
sigma <- 1.0
N <- 200
G <- 20   #number of groups
GroupStructure <- c(rep(p/G,G))
Ntest =1000
nrep = 500
PMSEstest1 = matrix(NA, nrow = nrep, ncol = 1)
PMSEstest2 = matrix(NA, nrow = nrep, ncol = 3)

colnames(PMSEstest2)<- c("5","10","20")

#CoData <- c(rep(1,100),rep(2,100),rep(3,100),rep(4,100),rep(5,100))
CoData <- c()
for (i in 1:G){
  CoData <- append(CoData,rep(i,p/G))
}
CoData

CoData <- factor(CoData, levels = unique(CoData))
library(ecpc)
Groupset <- createGroupset(CoData)
CoDataMatrix <- createZforGroupset(Groupset)
Z_all = list(CoDataMatrix)

g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,101] - 0.5)^2 + 10 * x[,102] + 10 * x[,3]
}
set.seed(1)
Xtest <- matrix(runif(Ntest*p),Ntest,p)
Ytest <- g(Xtest) + rnorm(Ntest,0,sigma)

for (j in 1:nrep){    # for loop representing simulated data sets
  print(paste("Sim","Iter",j,sep = " "))
  
  # simulate training data, either setting 1 or setting 2
  set.seed(j^3+239)
  
  X <- matrix(runif(N * p), N, p)
  Y <- g(X)+ rnorm(N, 0, sigma)
  res <- ecpc(Y=Y, X = X, Z = Z_all, Y2 = Ytest, X2 = Xtest, maxsel = c(5,10,20))
  PMSEstest1[j,1] = res$MSEecpc
  PMSEstest2[j,] = res$MSEPost
  remove(res); gc()
}
mean(PMSEstest1)
colMeans(PMSEstest2)
res = cbind.data.frame("noVarsel" = PMSEstest1,"Varsel"=PMSEstest2)


save(res,file = "ECPCresults_N=200,G=20.Rdata")

############### CORF ###############
####################################

## Set-up ##
p <- 500
sigma <- 1.0
N <- 200
G <- 5   #number of groups
GroupStructure <- c(rep(p/G,G))
Ntest =500
nrep = 500

g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,101] - 0.5)^2 + 10 * x[,102] + 10 * x[,3]
}
set.seed(1)
Xtest <- matrix(runif(Ntest*p),Ntest,p)
Ytest <- g(Xtest) + rnorm(Ntest,0,sigma)
DFtest <- data.frame(Ydf = Ytest, Xdf = Xtest)

remove(Xtest,Ytest); gc()


PMSEstest1 = matrix(NA, nrow = nrep, ncol = 1)
PMSEstest2 = matrix(NA, nrow = nrep, ncol = 1)


CoData <- c()
for (i in 1:G){
  CoData <- append(CoData,rep(i,p/G))
}


CoData <- factor(CoData, levels = unique(CoData))
CoData <- data.frame("Group" = CoData)
library(CoRF); library(randomForestSRC)


for (j in 1:nrep){    # for loop representing simulated data sets
  print(paste("Sim","Iter",j,sep = " "))
  
  # simulate training data, either setting 1 or setting 2
  set.seed(j^3+239)
  
  X <- matrix(runif(N * p), N, p)
  Y <- g(X)+ rnorm(N, 0, sigma)
  DFtrain <- data.frame(Ydf=Y,Xdf=X)
  
  set.seed(j)
  #Fit base RF
  baseRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",
                  importance=c("none"), nodesize = 5, setseed=1)
  ## Fitting co-data model ##
  VarsUsed <- baseRF$var.used
  #CoDataTrain$pvalsVUmc[is.na(CoDataTrain$pvalsVUmc)] <- mean(CoDataTrain$pvalsVUmc,na.rm=TRUE)
  
  # glm 
  CoDataModel  <- glm(VarsUsed/sum(VarsUsed) ~  Group, 
                      data=CoData,family='quasibinomial')
  
  # obtaining estimated sampling probabilities
  #predswt <- as.numeric(plogis(predict(CoDataModel)))
  predswt <- as.numeric(plogis(predict(CoDataModel)))
  predswt2 <- pmax(predswt-1/p,0)
  predswt2 <- predswt2/sum(predswt2)
  Mtry <- sum(predswt2!=0)/3
  print(c(predswt2[1],predswt[51],predswt2[101],predswt2[201]))
  ## Fit new RF with estimated sampling probabilities ##
  CoRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",
                xvar.wt = predswt2, mtry = Mtry,
                importance = c("none"), nodesize = 5, setseed = 1)
  preds1 <- predict.rfsrc(baseRF, newdata = DFtest, outcome = "train") # base RF
  preds2 <- predict.rfsrc(CoRF, newdata = DFtest, outcome = "train") # CoRF
  preds1 <- preds1$predicted # base RF
  preds2 <- preds2$predicted # CoRF
  PMSEstest1[j,1] = mean((preds1-DFtest$Ydf)^2) #base RF
  PMSEstest2[j,1] = mean((preds2-DFtest$Ydf)^2) #CoRF
  print(c(PMSEstest1[j,1],PMSEstest2[j,1]))
  remove(CoRF,baseRF,preds1,preds2); gc()
}

mean(PMSEstest1)
mean(PMSEstest2)
res = cbind.data.frame("baseRF" = PMSEstest1,"CoRF"=PMSEstest2)


save(res,file = "CoRFresults_N=200,G=5.Rdata")
