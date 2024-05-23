#### Try-outs ####
setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper")


p <- 500
sigma <- 1.0
N <- 100
G <- 5   #number of groups

CoDat = rep(1:G, rep(p/G,G))
CoDat = data.frame(factor(CoDat))
CoDat <- model.matrix(~0+., CoDat)
colnames(CoDat)  = paste0("Group ",1:G)

g <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,101] - 0.5)^2 + 10 * x[,102] + 10 * x[,3]
}

X <- matrix(runif(N * p), N, p)
#colnames(X)<-paste0("x",seq(1,p))
Y <- g(X)+ rnorm(N, 0, sigma)

source('EBcoBART_MainFunction.R')
load("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Application/DatForPaper.Rdata")
X<-dat$Xtrain
X<-as.matrix(X)

Y<-dat$Ytrain
CoDat<-dat$CoData
CoDat<-model.matrix(~0+.,CoDat)
X$IPI = factor(X$IPI)
dat=Dat_EBcoBART(X,CoData = CoDat)
X<-dat$X
CoDat<-dat$CoData

a=EBcoBART(Y=Y,X=X,CoData = CoDat, nIter = 15, model = "binary", 
           EB = T, Info = T, Seed = T,
           k = 2, alpha = .95, beta = 2)
a$SplittingProbs

