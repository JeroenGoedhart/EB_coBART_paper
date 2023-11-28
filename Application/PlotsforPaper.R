setwd("C:/Users/VNOB-0732/Desktop/R files/EB_coBART_paper/Application")

library(ggplot2); library(viridis); library(tidyr); library(gridExtra)

load("Results/PartialDependence1EBcoBART.Rdata")
load("Results/LearningCurveResults_EBcoBART_Performances_2_2_.95_Treeupdate.Rdata")
load("Results/LearningCurveResults_EBcoBART_EstHyperpars_2_2_.95_Treeupdate.Rdata")
LearnCurve_rigid <- results1
LearnCurve_flex <- results1
AUCs_def_rig
AUCs_EBco_rig
AUCs_def_rig <- LearnCurve_rigid$AUCs_def
AUCs_EBco_rig <- LearnCurve_rigid$AUCs_EBCo
AUCs_def_rig <- unlist(lapply(AUCs_def_rig, mean))
AUCs_EBco_rig <- unlist(lapply(AUCs_EBco_rig, mean))

AUCs_def_flex <- LearnCurve_flex$AUCs_def
AUCs_EBco_flex <- LearnCurve_flex$AUCs_EBCo
AUCs_def_flex <- unlist(lapply(AUCs_def_flex, mean))
AUCs_EBco_flex <- unlist(lapply(AUCs_EBco_flex, mean))

remove(results1,AUCs_EBco,AUCs_def)
Groupprobs <- results$Probs
WAIC <- results$waic
k <- results$k
alpha <- results$alpha

minima = apply(WAIC,1, function (x) which(x==min(x))) # WAIC minimum
minima
dim(Groupprobs[[1]])
EstProbs <- matrix(nrow=10, ncol = 140)
alpha_min <- c()
k_min <- c()

for (i in 1:10) {
  EstProbs[i,]<-Groupprobs[[i]][minima[i],]
  alpha_min[i] <- alpha[i,minima[i]]
  k_min[i] <- k[i,minima[i]]
  
}

dat <- cbind.data.frame("ntrain" = nseq, "Estimate" = k_min, "parameter" = "k")
dat1 <- cbind.data.frame("ntrain" = nseq, "Estimate" = alpha_min, "parameter" = "alpha")

dat2 <- rbind.data.frame(dat,dat1)

name <- paste("HypEstalpha","EBcoBart.pdf", sep = "_")
pdf(name, width=4,height=3)
p = ggplot(dat1, aes(x = ntrain, y = Estimate))  + coord_cartesian(ylim = c(0,max(dat1$Estimate))) +
  geom_point(size=1) +
  geom_line(linewidth=0.8) +
  theme_light() +
  scale_color_viridis(discrete = T, option = "G",begin = 0.2, end = 0.7) +
  theme(legend.title = element_blank(), legend.position = c(0.2,0.85)) + labs(x="Training set size",y="EB-estimate alpha")
p
dev.off()




alpha_min
k_min
Est_IPI = EstProbs[,140]


nmin = 20; nev = 10;
nseq <- round(seq(nmin,91,length.out=nev))
nseq
LearningCurve = cbind.data.frame("ntrain" = nseq,
                                 "EB-coBART" = AUCs_EBco_flex,
                                 "BART" = AUCs_def_flex)

LearningCurve1 = cbind.data.frame("ntrain" = nseq,
                                 "EBcoBART" = AUC_EBco,
                                 "IPI" = Est_IPI)
library(ggplot2); library(viridis)
name <- paste("LearningCurve","EBcoBart.pdf", sep = "_")
pdf(name, width=4,height=3)
p = ggplot(LearningCurve1, aes(x = ntrain, y = EBcoBART)) +
  geom_point(size=3) + coord_cartesian(ylim = c(0.5,0.75)) +
  geom_line(linewidth=0.8) +
  theme_light() +
  scale_color_viridis(discrete = T, option = "G",begin = 0.2, end = 0.7) +
  theme(legend.title = element_blank(), legend.position = c(0.2,0.85)) + labs(x="Training set size",y="AUC")
p

dev.off()
max_first  <- 0.75   # Specify max of first y axis
max_second <- 1 # Specify max of second y axis
min_first  <- 0.5   # Specify min of first y axis
min_second <- 0 # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale = (max_second - min_second)/(max_first - min_first)
shift = min_first - min_second

# Function to scale secondary axis
scale_function <- function(x){
  return (x*4-2)
}


# Function to scale secondary variable values
inv_scale_function <- function(x){
  return ((x+2)/4)
}
inv_scale_function(0)
a=scale_color_viridis(discrete = T, option = "G",begin = 0.2, end = 0.7)
a=viridis::mako(2,begin = 0.2, end = 0.7)
name <- paste("LearningVurveP_Perf_and_Weight","EBcoBart.pdf", sep = "_")
pdf(name, width=4,height=3)
pkpd <- ggplot(LearningCurve1, aes(x = ntrain, y = EBcoBART)) +
  geom_point(aes(color="Performance", shape="Performance"))+
  geom_line(aes(color="Performance"),linewidth=1) +
  geom_point(aes(y = inv_scale_function(IPI),color="Weight", shape = "Weight"))+
  geom_line(aes(y = inv_scale_function(IPI),color="Weight"),linewidth=1) +
  
  
  scale_y_continuous(limits = c(min_first, max_first), sec.axis = sec_axis(~scale_function(.), name="Weight of IPI")) +
  labs(x = "Training set size", y = "AUC") +
  theme_light() +
  scale_color_manual(name = "",values=a, labels = c("Performance","Weight")) +
  scale_shape_manual(name = "", values=c(17,19),labels = c("Performance","Weight"))+
  
  theme(legend.position = c(0.7,0.2)) 
print(pkpd)
dev.off()
scale_shape_m
names(LearningCurve)
lib
LearningCurve1 <- pivot_longer(LearningCurve, cols = c("EB-coBART","BART"), names_to = "Method", values_to = "Val")

PD$method[which(PD$method=="defaultBART")]<-"BART"
PD$method[which(PD$method=="EBcoBART")]<-"EB-coBART"

name <- paste("PartialDependence","EBcoBart.pdf", sep = "_")
pdf(name, width=4,height=3)
g = ggplot(PD, aes(x = IPI, y = average, color = method, shape = method)) +
  geom_errorbar(aes(ymin=average-sd, ymax=average+sd), width=.1) +
  geom_point(size = 3) +
  geom_line(linewidth=0.8) +
  theme_light() +
  scale_color_viridis(discrete = T, option = "G",begin = 0.2, end = 0.7)+
  theme(legend.title = element_blank(), legend.position = c(0.2,0.85)) + labs(x="IPI",y="Latent Response")
g

dev.off()




name <- paste("LearningCurve_EBcoBart.pdf", sep = "_")
pdf(name, width=4,height=3)
p = ggplot(LearningCurve1, aes(x = ntrain, y = Val, color = Method, shape = Method)) +
  geom_point(size=3) + coord_cartesian(ylim = c(0.5,0.75)) +
  geom_line(linewidth=0.8) +
  theme_light() +
  scale_color_viridis(discrete = T, option = "G",begin = 0.2, end = 0.7) +
  theme(legend.title = element_blank(),legend.position = c(0.2,0.85)) + labs(x="Training set size",y="AUC")
p

dev.off()
getwd()
#################################
gc()
load("Results/FinalModelEBCoBARTResults1_.95_2_2.Rdata")
library(ggplot2); library(tidyr); library(viridis); library(gridExtra)
waic <- results$WAICs
waic
Model <- results$CoData
alpha <- results$alpha
k <- results$k
alpha[9]
plot(alpha)
plot(k)
alpha <- results$hyperparameters

dat <- Model$data
preds <- predict(Model, type = "response")
preds
dat <- cbind.data.frame(dat,preds)

#dat$p.values = exp(-dat$p.values)/(1+exp(-dat$p.values))


CNV = dat[which(dat$Groups=="CopyNumber"),]
CNV1 = CNV[order(CNV$p.values, decreasing = F),]

Mut = dat[which(dat$Groups=="Mutation"),]
Mut1 = Mut[order(Mut$p.values,decreasing = F),]

Trans = dat[which(dat$Groups=="Translocation"),]
Trans1 = Trans[order(Trans$p.values,decreasing = F),]
remove(CNV,Mut,Trans)

CNV1$preds = cumsum(CNV1$preds)
Mut1$preds = cumsum(Mut1$preds)
Trans1$preds = cumsum(Trans1$preds)

dat1 = rbind.data.frame(CNV1,Mut1,Trans1,dat[140,])



waic

plot.waic = cbind.data.frame("Iteration"= seq(1,18,1),"WAIC" = waic)
library(ggplot2); library(viridis); library(cowplot)
cols = viridis(4)

name <- paste("EBcoBART_flex_True_GroupWeights.pdf")
pdf(name, width=4,height=3)
p1 <- ggplot(dat1, aes(x=p.values, y=preds, color=Groups,shape=Groups)) + # asking it to set the color by the variable "group" is what makes it draw three different lines
  geom_point(size=c(rep(2,139),3)) + theme_light() +
  #scale_color_viridis(discrete = T, option = "D") +
  labs(x="- logit(p-values)", y="Cumulative Weight", color = "Variable type") +
  theme(legend.position = c(0.2, 0.70)) +
  scale_colour_manual(name = "Group",
                      labels = c("CNV","Mutation","Translocation","IPI"),
                      values = cols) +   
  scale_shape_manual(name = "Group",
                     labels = c("CNV","Mutation","Translocation","IPI"),
                     values = c(15, 17, 18, 19))
p1  
dev.off()

name <- paste("EBcoBART_flex_True_WAIC.pdf")
pdf(name, width=4,height=3)
g1 <- ggplot(plot.waic, aes(x=Iteration, y=WAIC)) + # asking it to set the color by the variable "group" is what makes it draw three different lines
  geom_point() + theme_light() +
  scale_color_viridis(discrete = T, option = "A") +
  labs(x="Iteration", y="WAIC") +
  geom_vline(xintercept=which(waic==min(waic)),linetype=2,linewidth=1, col = "red")
g1
dev.off()

grid.arrange(p,g1,p1,g)
library(patchwork)
name <- paste("Figure3_Application.pdf")
pdf(name, width=8,height=6)
plot_grid(p,g1,p1,g,
          labels = c("a", "b", "c","d"))

dev.off()

name <- paste("Figure3_Application_Optional1.pdf")
pdf(name, width=8,height=3)
plot_grid(g1,p1,
          labels = c("a", "b"))

dev.off()

name <- paste("Figure3_Application_Optional2.pdf")
pdf(name, width=8,height=3)
plot_grid(p,g,
          labels = c("a", "b"))

dev.off()

name <- paste("Figure3_Application_Optional3.pdf")
pdf(name, width=4,height=3)
p

dev.off()



getwd()
