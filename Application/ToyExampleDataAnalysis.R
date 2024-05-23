#### Toy Example from caret package ####

setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Application")

LikelihoodCont <- function(Ypred, Y,sigma){
  # returns the bernoulli likelihoods for the sampled posterior parameters
  # output is a N x m matrix, with N the number of samples and m the number of mc samples
  
  loglik <- -(0.5*(1/sigma^2))*(sweep(Ypred,2,Y)^2)-.5*log(sigma^2)-.5*log(2*pi)
  return(loglik)
}

### loading relevant libraries ###
library(dbarts)
library(loo)
library(pROC)
#install.packages("caret")
library(caret)
library(multiridge)
data(BloodBrain)

#var(logBBB)
#View(cor(bbbDescr))
#1 performance plafond; 2 binary co-data significant covariates #2 bart gewichten als co-data
#4 spearman p-values

sds = apply(bbbDescr,MARGIN = 2,sd)
ids = which(sds<0.075)


bbbDescr=bbbDescr[,-ids] #only zeros




#### define co-data set ####
set.seed(434)
ids = sample(1:208,size=0.8*208,replace = F)

Codata=cbind.data.frame(bbbDescr[ids,],y=logBBB[ids])
sds = apply(Codata,MARGIN = 2,sd)
which(sds<0.05)



X = bbbDescr[-ids,]; Y = logBBB[-ids]
sds=apply(X,MARGIN = 2,sd)
which(sds<0.05)
remove(logBBB,bbbDescr)


#set.seed(1234)
#folds <- CVfolds(Y=Y, balance = F, kfold = 5, nrepeat = 3) ## folds for 5 times repeated 10 fold CV
save(X=X,Y=Y,CoDat = Codata,folds=folds,file = "dat_Toy.Rdata")

load("dat_Toy.Rdata")

### define Co-data for BART ###
# estimated var importances of BART model #
library(dbarts)
fit = bart(x.train = Codata[,-122], y.train = Codata$y, k = 2, power = 2, base = .95,
           sigdf = 10, sigest = 0.667*sd(Codata$y), sigquant = .75,
           ndpost = 20000L,                  # number of posterior samples
           nskip = 20000L,                   # number of "warmup" samples to discard
           nchain = 10L,                     # number of independent, parallel chains
           ntree = 50L,                      # number of trees per chain
           keeptrees = F,
           keepevery = 10L,
           seed = 1)

VarCount = fit$varcount
VarCount = colSums(fit$varcount)
VarProb = VarCount/sum(VarCount)

CoData_BART = data.frame(VarCount)

### Define codata for ecpc ###
library(glmnet)
cvglmnet = cv.glmnet(x=as.matrix(Codata[,-122]), y=Codata$y,alpha=0)
lam=cvglmnet$lambda.min
ridge = glmnet(x=as.matrix(Codata[,-122]), y=Codata$y,alpha=0,lambda = lam)
betas = ridge$beta
betas= matrix(betas)[,1]
CoData_ridge = cbind(rep(1,ncol(X)),betas)



### Define codata for CoRF ###
library(randomForestSRC)
baseRF <- rfsrc(y ~ .,data=Codata,ntree=2000,var.used="all.trees",importance=c("none"))
VarsUsed <- baseRF$var.used
CoData_rf = data.frame(VarsUsed)

save("BART" = CoData_BART,"Ridge" = CoData_ridge,"RF" =CoData_rf,file = "CodataMatrices_Toy.Rdata")



load("dat_Toy.Rdata"); load("CodataMatrices_Toy.Rdata");


### Fitting EB-coBART ###

nIter = 12; p = ncol(X)
# initialize model in rigid or flexible tree setting
model <- "rigid"
#model <- "flexible"
if (model=="rigid"){k = 1; base = .1; power = 4.0}
if (model=="flexible"){k=2; base = .95; power = 2.0}


# storage containers for results
GroupProbabilities <- vector("list", length = length(folds))
WAICs = matrix(NA, nrow = length(folds), ncol = nIter)


for (j in 1:length(folds)) {
  ids = folds[[j]]
  Xtrain = makeModelMatrixFromDataFrame(X[-ids,], drop = F)
  Ytrain <- Y[-ids]
  sigest = 0.667*sd(Ytrain)
  print(paste("Fold",j,sep = " "))
  EstimatedProbs <- matrix(NA, nrow = nIter+1, ncol = p)
  EstimatedProbs[1,] <- c(rep(1/p,p)) 
  row.names(EstimatedProbs) <- c(paste("Iter",0:(nIter), sep = " "))
  
  # storage containers for single fold
  WAICVector <- c()
  
  
  
  # initialization of EBcoBART for single fold
  probs <- c(rep(1/p,p))
  
  for (i in 1:nIter) {
    print(paste("iteration",i,sep = " "))
    fit <- bart(x.train = Xtrain, y.train = Ytrain,
                ndpost = 20000L,                  # number of posterior samples
                nskip = 20000L,                   # number of "warmup" samples to discard
                nchain = 10L,                     # number of independent, parallel chains
                ntree = 50L,                      # number of trees per chain
                keeptrees = F,
                keepevery = 10L,
                usequants = F,
                k = k, base = base, power = power, # hyperparameters tree
                sigest = sigest, sigdf = 10, sigquant = .75,
                splitprobs = probs,                # prob that variable is chosen for split
                combinechains = T,
                verbose = F,
                seed = 3*j^2+202+4*i) 
    
    # Codata Model
    VarsUsed <- colSums(fit$varcount)
    VarsUsed <- VarsUsed/sum(VarsUsed)
    
    
    coDataModel <- glm(VarsUsed ~ .,
                       data=CoData_BART,family=quasibinomial)
    
    probs <- predict(coDataModel, type="response", newdata = CoData_BART)
    probs <- unname(probs)
    probs[is.na(probs)] <- .0000000000000001
    EstimatedProbs[i+1,] <- probs
    remove(coDataModel,VarsUsed)
    
    # estimate waic
    Ypred = fit$yhat.train
    LogLikMatrix = LikelihoodCont(Ypred = Ypred, Y = Ytrain, sigma = fit$sigma)
    WAICVector[i] <- suppressWarnings(waic(LogLikMatrix)$estimates[3,1])
    remove(LogLikMatrix,Ypred,fit)
    
    
  }
  
  WAICs[j,] <- WAICVector
  GroupProbabilities[[j]] <- EstimatedProbs
  remove(Xtrain,Ytrain)
  

}
plot(WAICs[12,])
View(GroupProbabilities[[1]])

results=list(Estprobs = GroupProbabilities,WAIC = WAICs, CVfolds = folds)
getwd()
save(results, file="EBcoBARTfit_rigid_ToyExample.Rdata")

load("EBcoBARTfit_rigid_ToyExample.Rdata")

WAICs <- results$WAIC
EstProbs <- results$Estprobs


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

maxTrain = apply(WAICs,1, function (x) which.min(x)) # find minimum WAICS for each fold
maxTrain
maxTrain = apply(WAICs,1, function (x) minimum(x))
maxTrain

# determine estimated covariate weights at minimum WAIC
MinProbs = matrix(NA, nrow = length(maxTrain), ncol = ncol(X))
for (i in 1:length(maxTrain)){
  MinProbs[i,] <- EstProbs[[i]][maxTrain[i],]
}
remove(EstProbs,WAICs,i,j)
#model <- "flexible"
model <- "rigid"
if (model=="rigid"){k = 1; base = .1; power = 4.0}
if (model=="flexible"){k=2; base = .95; power = 2.0}



MSE_coBART_rig <- c()
MSE_defBART_rig <- c()

gc()

library(dbarts)



for (j in 1:length(folds)) {
  
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = makeModelMatrixFromDataFrame(X[-ids,], drop = F)
  Ytrain <- Y[-ids]; print(var(Ytrain))
  sigest = 0.667*sd(Ytrain)
  Xtest = makeModelMatrixFromDataFrame(X[ids,], drop = F)
  Ytest <- Y[ids]
  
  fit <- bart(x.train = Xtrain, y.train = Ytrain,
              x.test = Xtest,
              ndpost = 20000L,   # number of posterior samples
              nskip = 20000L, # number of "warmup" samples to discard
              nchain = 10L,   # number of independent, parallel chains
              ntree = 50L,    # number of trees per chain
              keeptrees = F,
              verbose = F,
              usequants = F,
              k = k, base = base, power = power, # hyperparameters tree
              sigest = sigest, sigdf = 10, sigquant = .75,
              splitprobs = MinProbs[j,], keepsampler = F) 
  
  ypred <- fit$yhat.test
  ypred <- colMeans(ypred)
  MSE_coBART_rig[j] <- mean((ypred - Ytest)^2)
  
  
  remove(fit,ypred)
  
  fit1 <- bart(x.train = Xtrain, y.train = Ytrain,
               x.test = Xtest,
               ndpost = 20000L,   # number of posterior samples
               nskip = 20000L, # number of "warmup" samples to discard
               nchain = 10L,   # number of independent, parallel chains
               ntree = 50L,    # number of trees per chain
               keeptrees = F,
               verbose = F,
               usequants = F,
               k = k, base = base, power = power, # hyperparameters tree
               sigest = sigest, sigdf = 10, sigquant = .75,
               splitprobs = c(), keepsampler = F) 
  
  ypred1 <- fit1$yhat.test
  ypred1 <- colMeans(ypred1)
  MSE_defBART_rig[j] <- mean((ypred1 - Ytest)^2)
  remove(fit1,ypred1)
  remove(Ytrain,Ytest,Xtrain,Xtest)
}

mean(MSE_coBART_rig);mean(MSE_defBART_rig)
perfres = cbind.data.frame(perfres,'EBcoRig' = MSE_coBART_rig, 'BARTRig' = MSE_defBART_rig)
colMeans(perfres)
save(perfres,file="Perf_ToyExample.Rdata")
load("Perf_ToyExample.Rdata")                      



### Fit ecpc ###
library(ecpc)
MSE_ridge <- c()
MSE_ecpc <- c()
Yte <- c()
for (j in 1:length(folds)) {
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = X[-ids,]
  Ytrain <- Y[-ids]
  
  Xtest = X[ids,]
  Ytest <- Y[ids]
  Yte[j]<- var(Ytest)
  fit <- ecpc(Y = Ytrain, X = as.matrix(Xtrain), Z = list(CoData_ridge), postselection = F,
              X2 = as.matrix(Xtest), Y2 = Ytest, maxsel = c(5,10,50,75))
  fit$MSEecpc
  fit$MSEPost
  print(fit$gamma)
  MSE_ridge[j] = fit$MSEridge
  MSE_ecpc[j] = fit$MSEecpc
  remove(Xtrain,Ytrain,Xtest,Ytest,fit)
}
mean(MSE_ridge); mean(MSE_ecpc)
mean(Yte)
load("Perf_ToyExample.Rdata")
perfres = cbind.data.frame(perfres,"Ridge" = MSE_ridge, "ECPC" = MSE_ecpc,"VarY" = Yte)
save(perfres,file="Perf_ToyExample.Rdata")


### Fit CoRF ###
library(randomForestSRC)
MSE_rf <- c()
MSE_CorF <- c()
for (j in 1:length(folds)) {
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = X[-ids,]; Ytrain <- Y[-ids]
  Xtest = X[ids,]; Ytest <- Y[ids]
  DFtrain <- data.frame(Ydf=Ytrain,Xdf=Xtrain)
  DFtest <- data.frame(Ydf = Ytest, Xdf = Xtest)
  ### Fitting baseRF ###
  baseRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",
                  importance=c("none"),splitrule = 'mse')
  preds1 <- predict.rfsrc(baseRF, newdata = DFtest)$predicted # base RF
  MSE_rf[j]<- mean((preds1 - Ytest)^2)
  
  
  ## Fitting co-data model ##
  VarsUsed <- baseRF$var.used
  
  CoDataModel <-  glm(VarsUsed/sum(VarsUsed) ~ 1 + VarsUsed, 
                      data=CoData_rf,family='quasibinomial')
  predswt <- as.numeric(plogis(predict(CoDataModel)))
  P <- length(predswt)
  predswt2 <- pmax(predswt-1/P,0)
  Mtry <- ceiling(sqrt(sum(predswt2!=0)))
  ### Fitting Corf ###
  CoRF <- rfsrc(Ydf ~ .,data=DFtrain,ntree=2000,var.used="all.trees",
                xvar.wt = predswt2, mtry = Mtry,splitrule = "mse")
  
  ## Estimating performance for base RF and CoRF ##
  preds2 <- predict.rfsrc(CoRF, newdata = DFtest, outcome = "train")$predicted # CoRF
  MSE_CorF[j]<- mean((preds2 - Ytest)^2)
  
  remove(preds1,preds2,baseRF,CoRF,Xtrain,Ytrain,Xtest,Ytest,
         CoDataModel,predswt,predswt2,Mtry,VarsUsed,ids,DFtest,DFtrain)
}
mean(MSE_CorF);mean(MSE_rf)
load("Perf_ToyExample.Rdata")
perfres = cbind.data.frame(perfres,"RF" = MSE_rf, "CoRF" = MSE_CorF)
save(perfres,file="Perf_ToyExample.Rdata")
colMeans(perfres)

Rsquar = apply(perfres,2,function (x) 1-x/perfres$VarY)
colMeans(Rsquar)
colMeans(perfres)

warnings()



perfres = cbind.dataframe('EBcoFlex' = MSE_coBART, 'BARTflex' = MSE_defBART,
                          'Ridge' = MSE_ridge,'ECPC' = MSE_ecpc,
                          'RF' = MSE_rf, 'CoRF' = MSE_CorF)



mean(MSE_coBART);mean(MSE_defBART);mean(MSE_ridge);mean(MSE_ecpc);mean(MSE_rf); mean(MSE_CorF);mean(Yte)












### estimate perf on whole data set ###
set.seed(1234)
folds <- CVfolds(Y=logBBB, balance = F, kfold = 10, nrepeat = 3) ## folds for 5 times repeated 10 fold CV
model <- "flexible"
probs <- c(rep(1/ncol(bbbDescr),ncol(bbbDescr)))
if (model=="rigid"){k = 1; base = .1; power = 4.0}
if (model=="flexible"){k=2; base = .95; power = 2.0}
MSE_BART <- c()

for (j in 1:length(folds)) {
  
  print(paste("Fold",j,sep = " "))
  ids = folds[[j]]
  Xtrain = makeModelMatrixFromDataFrame(bbbDescr[-ids,], drop = F)
  Ytrain <- logBBB[-ids]; print(var(Ytrain))
  sigest = 0.667*sd(Ytrain)
  Xtest = makeModelMatrixFromDataFrame(bbbDescr[ids,], drop = F)
  Ytest <- logBBB[ids]
  
  fit <- bart(x.train = Xtrain, y.train = Ytrain,
              x.test = Xtest,
              ndpost = 10000L,   # number of posterior samples
              nskip = 30000L, # number of "warmup" samples to discard
              nchain = 10L,   # number of independent, parallel chains
              ntree = 50L,    # number of trees per chain
              keeptrees = F,
              verbose = F,
              usequants = F,
              k = k, base = base, power = power, # hyperparameters tree
              sigest = sigest, sigdf = 10, sigquant = .75,
              splitprobs = probs, keepsampler = F) 
  
  ypred <- fit$yhat.test
  ypred <- colMeans(ypred)
  MSE_BART[j] <- mean((ypred - Ytest)^2)
}
mean(MSE_BART)



### define Co-data: estimate p-values from Co-data set ###
varnames=names(X)
xnam <- paste0("x", 1:length(varnames))
xnam = append(xnam,"y")
names(Codata) <- xnam
p.values_spear <- c()
p.values_pears <- c()
for (i in 1:length(varnames)){
  fmla <- as.formula(paste("y ~ ", paste(names(Codata[i]))))
  a=lm(formula = fmla, data = Codata)
  print(coef(summary(a)))
  p.value = coef(summary(a))[-1,4] 
  p.values_pears[i] = p.value
  b=cor.test(x=Codata[,i],y=Codata$y,method = "spearman")
  p.values_spear[i] = b$p.value
  
  #gc()
}

View(cbind(p.values_pears,p.values_spear))


p.values.adjust <- p.adjust(p.values_spear, method = "BH", n = length(p.values_spear))

which(p.values.adjust<0.05)

CoData = data.frame(pval=-log(p.values.adjust/(1-p.values.adjust)))
#CoData=data.frame(EffEst = p.values)
rownames(CoData)=varnames;remove(Codata)

remove(fmla,i,ids,p.values,p.values.adjust,xnam)



