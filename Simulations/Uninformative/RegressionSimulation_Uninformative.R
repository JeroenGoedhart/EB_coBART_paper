gc()
setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Simulations/Uninformative")


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
Ntest = 500
G <- 5   #number of groups
GroupStructure <- c(rep(p/G,G))


g <- function(x){
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 10 * x[,5]
}


dat = vector(type='list',length = G)
for (i in 1:G){
  dat[[i]] = matrix(runif(N * GroupStructure[i]), N, GroupStructure[i])
}
y=lapply(dat, function (x) g(x) + rnorm(N,0,sigma))
X = do.call(cbind,dat); remove(dat)
Y = Reduce("+",y)


# G = 5
g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 10 * x[,4] + 10 * x[,5] +
    10 * sin(pi * x[,101] * x[,102]) + 20 * (x[,103] - 0.5)^2 + 10 * x[,104] + 10 * x[,105] +
    10 * sin(pi * x[,201] * x[,202]) + 20 * (x[,203] - 0.5)^2 + 10 * x[,204] + 10 * x[,205] +
    10 * sin(pi * x[,301] * x[,302]) + 20 * (x[,303] - 0.5)^2 + 10 * x[,304] + 10 * x[,305] +
    10 * sin(pi * x[,401] * x[,402]) + 20 * (x[,403] - 0.5)^2 + 10 * x[,404] + 10 * x[,405]
}

# G = 10
g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 5 * x[,4] 
  10 * sin(pi * x[,51] * x[,52]) + 20 * (x[,53] - 0.5)^2 + 5 * x[,54] 
  10 * sin(pi * x[,101] * x[,102]) + 20 * (x[,103] - 0.5)^2 + 5 * x[,104] 
  10 * sin(pi * x[,151] * x[,152]) + 20 * (x[,153] - 0.5)^2 + 5 * x[,154] 
  10 * sin(pi * x[,401] * x[,402]) + 20 * (x[,403] - 0.5)^2 + 5 * x[,404] 
}


# G = 15

# G = 20

# G = 50





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

name <- paste("Uninformative",G,N, k, base, power,"EBBart_WAIC.Rdata", sep = "_")
results <- list(GroupProbs = GroupProbabilities,TestPerformance = PMSEstest, waic = WAICs)

save(results, file = name)
