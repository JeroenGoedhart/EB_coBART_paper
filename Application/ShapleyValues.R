setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Application")

load("Results/EBCoBARTfit_flexible_TRUE_.Rdata")
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
library(shapviz)
library(fastshap)
library(viridis)
library(ggplot2)
library(cowplot)
### EBcoBART
fit <- bart(x.train = Xtrain, y.train = Ytrain,
            #x.test = Xtest,
            ndpost = 5000L,   # number of posterior samples
            nskip = 12000L, # number of "warmup" samples to discard
            nchain = 5L,   # number of independent, parallel chains
            ntree = 50L,    # number of trees per chain
            keeptrees = TRUE,
            verbose = F,
            usequants = F,
            k = k_min, base = base_min, power = 2, # hyperparameters tree
            splitprobs = probsmax, keepsampler = T, seed = 22333) 

##### SHapley values

pfun <- function(object, newdata) {
  colMeans(predict(object, newdata = newdata))
}
set.seed(42)
shapall <-
  explain(fit,
          X = Xtrain,
          pred_wrapper = pfun,
          nsim = 5,
          newdata = Xtrain,
          adjust = T, shap_only = F)



sv <- shapviz(shapall)

### standard BART
fit <- bart(x.train = Xtrain, y.train = Ytrain,
            #x.test = Xtest,
            ndpost = 5000L,   # number of posterior samples
            nskip = 12000L, # number of "warmup" samples to discard
            nchain = 5L,   # number of independent, parallel chains
            ntree = 50L,    # number of trees per chain
            keeptrees = TRUE,
            verbose = F,
            usequants = F,
            k = 2, base = .95, power = 2, # hyperparameters tree
            splitprobs = c(), keepsampler = T, seed = 22333) 

##### SHapley values

pfun <- function(object, newdata) {
  colMeans(predict(object, newdata = newdata))
}
set.seed(42)
shapall1 <-
  explain(fit,
          X = Xtrain,
          pred_wrapper = pfun,
          nsim = 5,
          newdata = Xtrain,
          adjust = T, shap_only = F)



sv1 <- shapviz(shapall1)

col = viridis(2, option = "F",begin = 0.1, end = 0.5)
#### plotting
shapEBco <- colMeans(abs(shapall$shapley_values))
shapBART <- colMeans(abs(shapall1$shapley_values))
shapEBco <-sort(shapEBco,decreasing = T)
shapBART <- sort(shapBART,decreasing = T)

ids = which(names(shapEBco)[1:10] %in% names(shapBART)[1:10])
ids
shapEBco <- shapEBco[ids] 
ids = which(names(shapBART)[1:10] %in% names(shapEBco))
shapBART <- shapBART[ids]
shapEBco
shapBART<- shapBART[names(shapEBco)]
res=cbind.data.frame(shapBART,shapEBco,names(shapEBco))
colnames(res)<-c("BART","EB-coBART","Variable")
res1 <- pivot_longer(res, cols = c("BART","EB-coBART"), names_to = "Method", values_to = "Val")
res1
res1$Method <- factor(res1$Method, levels = c("BART","EB-coBART"))
res1$Variable<-factor(res1$Variable)
name <- paste("Figure4_Application.pdf")
pdf(name, width=4,height=3)
p<-ggplot(data=res1, aes(x=Variable, y=Val, fill = Method,density=Method ,angle=Method)) +
  geom_bar(stat="identity", color="black", position=position_dodge2(reverse=T)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "F",begin = 0.1, end = 0.5) +
  
  theme(legend.title = element_blank(), legend.position = c(0.8,0.2)) + labs(y="Average absolute Shapley value")

p + coord_flip()
dev.off()
?geom_bar()
p <- sv_importance(sv, kind = "beeswarm", max_display = 5,
              fill = col[2], bar_width = 0.2,
              viridis_args=list(begin = 0.25, end = 0.85,option = "E"))
p
p<- p+theme_light()
p
q <- sv_importance(sv1, kind = "bar", max_display = 5,
                   fill = col[1], bar_width = 0.2,
                   viridis_args=list(begin = 0.25, end = 0.85,option = "E"))
q<- q + theme_light()

library(patchwork)
name <- paste("Figure4_Application.pdf")
pdf(name, width=8,height=3)
plot_grid(p,q,
          labels = c("EB-coBART", "BART"))

dev.off()
getwd()
