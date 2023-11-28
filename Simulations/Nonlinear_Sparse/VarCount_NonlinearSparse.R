setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Simulations/Nonlinear_Sparse")
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
library(dbarts)
library(loo)

p <- 500
sigma <- 1.0
N <- 100
G <- 5   #number of groups
GroupStructure <- c(rep(p/G,G))

g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,101] - 0.5)^2 + 10 * x[,102] + 10 * x[,3]
}

#model <- "rigid"
model <- "flexible"



if (model=="rigid"){k = 1; base = .1; power = 4.0}
if (model=="flexible"){k=2; base = .95; power = 2.0}
nu <- 10 ; quant <- .75 # Error variance hyperparameters



### compute est. probs for EBcoBART
load("Results/Friedman_5_100_2_0.95_2_EBBart_WAIC_Final.Rdata")
GroupProbs <- results$GroupProbs
waic <- results$WAIC
remove(results); gc()
minima = apply(waic,1, function (x) minimum(x)) # WAIC minimum
EstProbs <- matrix(nrow=500, ncol = G)

for (i in 1:500) {
  EstProbs[i,]<-GroupProbs[[i]][minima[i],]
  
}

nrep = 500
VarCountsEBco = vector("list", length = nrep)
VarCountsDef = vector("list", length = nrep)

for (j in 1:nrep){    # for loop representing simulated data sets
  print(paste("Sim","Iter",j,sep = " "))
  
  # simulate training data, either setting 1 or setting 2
  set.seed(j^3+239)
  X <- matrix(runif(N * p), N, p)
  colnames(X) <- paste0("x",seq(1,500))
  Y <- g(X)+ rnorm(N, 0, sigma)
  
  # hyperparameter initialization
  sigest <- sd(Y)*0.667
  fitDefault <- bart(x.train = X, y.train = Y, # training data
              ndpost = 10000L,   # number of posterior samples
              nskip = 5000L, # number of "warmup" samples to discard
              nchain = 5L,   # number of independent, parallel chains
              ntree = 50L,    # number of trees per chain
              keeptrees = F,
              verbose = F,
              k = k, base = base, power = power, # hyperparameters tree
              sigest = sigest,sigdf = nu, sigquant = quant,  # hyperparameters error variance
              splitprobs = c()   # hyperparameter that will be updated using EB and co-data
  )
  VarCountsDef[[j]] <- sort(colSums(fitDefault$varcount),decreasing = T)
  
  
  weights <- c()
  
  for (i in 1:G) {
    weights <- append(weights,rep(EstProbs[j,i],p/G))
    
  }
  weights
  
  fitEBco <- bart(x.train = X, y.train = Y, # training data
                     ndpost = 10000L,   # number of posterior samples
                     nskip = 5000L, # number of "warmup" samples to discard
                     nchain = 5L,   # number of independent, parallel chains
                     ntree = 50L,    # number of trees per chain
                     keeptrees = F,
                     verbose = F,
                     k = k, base = base, power = power, # hyperparameters tree
                     sigest = sigest,sigdf = nu, sigquant = quant,  # hyperparameters error variance
                     splitprobs = weights   # hyperparameter that will be updated using EB and co-data
  ) 
  VarCountsEBco[[j]] <- sort(colSums(fitEBco$varcount),decreasing = T)
}
VarCountsDef[[2]][1:20]
VarCountsEBco[[2]][1:20]
results = list("DefBart" = VarCountsDef, "EBcoBART" = VarCountsEBco)
name <- paste("VarSel_Friedman",model,G,N,".Rdata", sep = "_")

save(results,file = name)

