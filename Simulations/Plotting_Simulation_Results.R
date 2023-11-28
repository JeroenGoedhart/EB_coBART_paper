setwd("C:/Users/VNOB-0732/Desktop/R files/Bayesian Additive Regression Trees/Grouping Information/GroupingBART/SimulationsEBBart")

gc()
### settings
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



model <- "flexible"
model <- "rigid"
Setting <- "G = 5, N = 100"

### results
PMSEstest = results$TestPerf
WAIC = results$WAIC
GroupProbabilities = results$GroupProbs

remove(results); gc()

#PlotResults <- c()

minima = apply(WAIC,1, function (x) minimum(x)) # WAIC minimum





EbcoBART <- c()
nrep <- length(minima)
for (i in 1:nrep){
  EbcoBART[i] = PMSEstest[i,minima[i]]
  #RatioPMSEtest1[i] = PMSEstest[i,minima1[i]]/PMSEstest[i,minima[i]]
}
mean(PMSEstest[,1]); mean(EbcoBART)


RatioPMSEtest <- c()
nrep <- length(minima)
for (i in 1:nrep){
  RatioPMSEtest[i] = PMSEstest[i,minima[i]]/PMSEstest[i,1]
  #RatioPMSEtest1[i] = PMSEstest[i,minima1[i]]/PMSEstest[i,minima[i]]
}
length(which(RatioPMSEtest>1))
median(RatioPMSEtest)
res <- cbind.data.frame("Perf" = RatioPMSEtest, "Model" = c(rep(model,500)), "Setting" = Setting)

#PlotResults <- c()
PlotResults <- rbind(PlotResults,res)


PlotResults$Model <- factor(PlotResults$Model, levels = c(unique(PlotResults$Model)))
PlotResults$Setting <- factor(PlotResults$Setting, levels = unique(PlotResults$Setting))
save(PlotResults, file = "PlotResults_PMSERatio_uninformative.Rdata")

PlotResults$Setting[501:1000]=unique(PlotResults$Setting)[2]
PlotResults$Setting[1501:2000]=unique(PlotResults$Setting)[1]

setwd("C:/Users/VNOB-0732/Desktop/R files/Bayesian Additive Regression Trees/Grouping Information/GroupingBART/SimulationsEBBart/Final Results")

load("PlotResults_PMSERatio.Rdata")
library(ggplot2); library(viridis)

name <- paste("DecreasePerf_UninformativeFriedman","EBBart_WAIC.pdf", sep = "_")
pdf(name, width=5,height=2.5)
bp <- ggplot(PlotResults, aes(x=Model, y=Perf, group=Model))  + coord_cartesian(ylim =  c(0.9, 1.1)) +
  geom_boxplot(aes(fill=Model),fatten =1 ,outlier.shape = NA) +
  #scale_y_continuous(breaks=seq(0.8,2.1,0.2)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "G",begin = 0.4, end = 1)+
  theme(legend.title = element_blank()) + labs(x="",y="PMSE Ratio") +
  
  geom_hline(yintercept=1,linetype=2,linewidth=1)
bp
#bp + facet_grid(. ~ Setting)
dev.off()
getwd()
#
#########################################
#PlotResults1 <- c()
res1 = cbind.data.frame("Avg. PMSE" = colMeans(PMSEstest),"Iteration" = c(seq(1,8)),"Model" = c(rep(model,8)), "Setting" = c(rep(Setting,8)))
res1$Model = factor(res1$Model, levels = unique(res1$Model))
res1$Setting = factor(res1$Setting, levels = unique(res1$Setting))
res1 <- unname(res1)
PlotResults1 <- unname(PlotResults1)
name = c("Avg. PMSE"."Iteration","Model","Setting")
names(PlotResults1) <- name
names(res1) <- name

PlotResults1 <- rbind(PlotResults1,res1)

PlotResults1$Model <- factor(PlotResults1$Model, levels = c(unique(PlotResults$Model)))
PlotResults1$Setting <- factor(PlotResults1$Setting, levels = unique(PlotResults$Setting))
names(PlotResults1)<- c("PMSE","Iteration","Model","Setting")
save(PlotResults1, file = "PlotResults_PMSE_Iter.Rdata")

setwd("C:/Users/VNOB-0732/Desktop/R files/Bayesian Additive Regression Trees/Grouping Information/GroupingBART/SimulationsEBBart/Final Results")

load("PlotResults_PMSE_Iter.Rdata")

library(viridis)
cols = viridis(2, direction = 1, option = "H")
name <- paste("PerfvsIter_Friedman","EBBart_WAIC.pdf", sep = "_")
pdf(name, width=7.5,height=2)
bp <- ggplot(PlotResults1, aes(x=Iteration, y= PMSE, group=Model, color=Model)) +
  geom_line(size=1) + geom_point() + theme_light() +  #coord_cartesian(ylim =  c(-0.2, 0.6)) +
  scale_color_viridis(discrete = T, option = "G",begin = 0.4, end = 1)+
  theme(legend.title = element_blank()) + labs(x="Iteration",y="Avg. PMSE")
bp + facet_grid(. ~ Setting)
dev.off()


############################################



#res_varsel <- c()
#res_varsel <- c()
G = 20

model <- "flexible"
model <- "rigid"
Setting <- "G = 5, N = 200"

### results
WAIC = results$WAIC
GroupProbabilities = results$GroupProbs
minima = apply(WAIC,1, function (x) minimum(x)) # WAIC minimum
EstProbs <- matrix(nrow=500, ncol = G)

for (i in 1:500) {
  EstProbs[i,]<-GroupProbabilities[[i]][minima[i],]
  
}

names(EstProbs) = paste0("G", seq(1,G,1),sep = "")
EstProbs <- cbind.data.frame(EstProbs,"Model"=rep(model,500),"Setting" = rep(Setting,500))
res_varsel <- rbind(res_varsel,EstProbs)


res_varsel$Model<- factor(res_varsel$Model, levels = unique(res_varsel$Model))
res_varsel$Setting<- factor(res_varsel$Setting, levels = unique(res_varsel$Setting))

save(res_varsel,file = "results_VarSel_Friedman_G=20_Friedman.Rdata")
getwd()
load("results_VarSel_Friedman_G=20.Rdata")
G=20
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
name <- c()
for (i in 1:G) {
  name[i] <- as.character(i)
  
}
name <- append(name,c("Model","Setting"))

names(res_varsel) <- name

res_varsel1 <- pivot_longer(res_varsel, cols = c(1,2,3,4,5,6,7,8,9,10), names_to = "Group", values_to = "Val")
res_varsel1 <- pivot_longer(res_varsel, cols = c(1,2,3,4,5), names_to = "Group", values_to = "Val")

res_varsel1$Group <- factor(res_varsel1$Group, levels = unique(res_varsel1$Group))
library(viridis);library(ggplot2)

getwd()
name <- paste("GroupWeights_Friedman_G20_EBBart_WAIC.pdf", sep = "_")
pdf(name, width=7,height=3)
bp <- ggplot(res_varsel1, aes(x=Group, y=Val)) +
  geom_boxplot(aes(fill=Model),fatten =1, outlier.shape = NA) + #coord_cartesian(ylim =  c(0.9, 1.3)) +
  #scale_y_continuous(breaks=seq(0.8,2.1,0.2)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "G",begin = 0.4, end = 1)+
  theme(legend.title = element_blank(), legend.position = c(0.92,0.82)) + labs(x="Group",y="Group Weights") +
  geom_hline(yintercept=0.2,linetype=3,linewidth=1)
bp
bp1 = bp + facet_grid(. ~ Setting)
bp1
gp <- ggplot(res_varsel1, aes(x=Group, y=Val)) +
  geom_boxplot(aes(fill=Model),fatten =1, outlier.shape = NA) + #coord_cartesian(ylim =  c(0.9, 1.3)) +
  #scale_y_continuous(breaks=seq(0.8,2.1,0.2)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "G",begin = 0.4, end = 1)+
  theme(legend.title = element_blank(), legend.position = c(0.92,0.82)) + labs(x="Group",y="Group Weights") +
  geom_hline(yintercept=0.05,linetype=3,linewidth=1)
gp
gp1 = gp + facet_grid(. ~ Setting)
gp1





library(cowplot)
name <- paste("Figure1_Simulations.pdf", sep = "_")
pdf(name, width=8,height=6)
plot_grid(bp1, gp1, ncol = 1,
          labels = c("a", "b"))
dev.off()


############################################################33
#### linear ####

#################################################################################

PlotResults_LinDense <- c()
PlotResults1_LinDense <- c()


model <- "flexible"
model <- "rigid"
Setting <- "N = 200"
PMSEstest = results$TestPerformance
WAIC = results$WAIC
remove(results); gc()


minima = apply(WAIC,1, function (x) which(x == min(x))) # WAIC minimum
RatioPMSEtest <- c()

nrep <- length(minima)
for (i in 1:nrep){
  RatioPMSEtest[i] = PMSEstest[i,1]/PMSEstest[i,minima[i]]
  #RatioPMSEtest1[i] = PMSEstest[i,minima1[i]]/PMSEstest[i,minima[i]]
}


res <- cbind.data.frame("Perf" = RatioPMSEtest, "Model" = c(rep(model,500)), "Setting" = Setting)
PlotResults_LinDense <- rbind(PlotResults_LinDense,res)

res1 = cbind.data.frame("Avg. PMSE" = colMeans(PMSEstest),"Iteration" = c(seq(1,10)),"Model" = c(rep(model,10)), "Setting" = c(rep(Setting,10)))
PlotResults1_LinDense <- rbind(PlotResults1_LinDense,res1)

PlotResults1_LinDense$Model<- factor(PlotResults1_LinDense$Model, levels = c(unique(PlotResults1_LinDense$Model)))
PlotResults1_LinDense$Setting <- factor(PlotResults1_LinDense$Setting, levels = unique(PlotResults1_LinDense$Setting))
names(PlotResults1_LinDense)<- c("PMSE","Iteration","Model","Setting")
PlotResults_LinDense$Perf = 1/PlotResults_LinDense$Perf
save(PlotResults1_LinDense, file = "PlotResults_PMSE_Iter_LinDense.Rdata")
save(PlotResults_LinDense, file = "PlotResults_PMSERatio_LinDense.Rdata")

load("PlotResults_PMSERatio_LinDense.Rdata")
load("PlotResults_PMSE_Iter_LinDense.Rdata")

library(ggplot2)
7.5/2
name <- paste("RatioPerf_LinDense","EBBart_WAIC.pdf", sep = "_")
pdf(name, width=3.75,height=2)
bp <- ggplot(PlotResults_LinDense, aes(x=Model, y=Perf, group=Model)) +
  geom_boxplot(aes(fill=Model),outlier.shape = NA) + coord_cartesian(ylim =  c(0.75, 1.07)) +
  #scale_y_continuous(breaks=seq(0.8,2.1,0.2)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "B",begin = 0.2, end = 0.7)+
  theme(legend.title = element_blank()) + labs(x="",y="PMSE Ratio") +
  geom_hline(yintercept=1,linetype=2,size=1)
bp + facet_grid(. ~ Setting)
dev.off()


name <- paste("PerfvsIter_LinDense","EBBart_WAIC.pdf", sep = "_")
pdf(name, width=3.75,height=2)
bp <- ggplot(PlotResults1_LinDense, aes(x=Iteration, y= PMSE, group=Model, color=Model)) +
  geom_line(size=1) + geom_point() + theme_light() +
  scale_color_viridis(discrete = T, option = "B",begin = 0.2, end = 0.7)+
  theme(legend.title = element_blank()) + labs(x="Iteration",y="Avg. PMSE")
bp + facet_grid(. ~ Setting)
dev.off()


#res_varsel <- c()

GroupProbabilities = results$GroupProbs
WAIC = results$WAIC
minima = apply(WAIC,1, function (x) which(x == min(x))) # WAIC minimum
model <- "flexible"
model <- "rigid"
Setting <- "N = 200"
EstProbs <- matrix(nrow=500, ncol = 500)

for (i in 1:500) {
  EstProbs[i,]<-GroupProbabilities[[i]][minima[i],]
  
}
remove(WAIC,GroupProbabilities,results)
EstProbs <- data.frame(EstProbs)
names(EstProbs) = paste0(seq(1,500,1),sep = "")
EstProbs <- cbind.data.frame(EstProbs,"Model"=rep(model,500),"Setting" = rep(Setting,500))

res_varsel <- rbind(res_varsel,EstProbs)

res_varsel$Model<- factor(res_varsel$Model, levels = unique(res_varsel$Model))
res_varsel$Setting<- factor(res_varsel$Setting, levels = unique(res_varsel$Setting))

save(res_varsel,file = "results_VarSel_LinDense.Rdata")
gc()
setwd("C:/Users/VNOB-0732/Desktop/R files/Bayesian Additive Regression Trees/Grouping Information/GroupingBART/SimulationsEBBart/Final Results")
load("results_VarSel_LinDense.Rdata")
nrow(res_varsel)
res_varsel <- res_varsel[,1:500]

G = 10

ids <- matrix(ncol=2,nrow = G)
for (i in 1:G){
  ids[i,1]<- (1+50*(i-1))
  ids[i,2]<-50*i 
}


res_tot <- c()
for (i in 1:G) {
  start = ids[i,1]
  end = ids[i,2]
  res = rowMeans(res_varsel[,start:end])
  res_tot <- cbind(res_tot,res)
}
res_tot <- cbind.data.frame(res_tot,res_varsel[,501:502])  
names(res_tot)[1:10] <- c("1","2","3","4","5","6","7","8","9","10")



res_tot 




res_tot <- c()
for (i in 1:nrow(res_varsel)) {
  res <- c()
  for (j in 1:G) {
    cols = :(100*j)
    res[j] <- mean(res_varsel[i,cols])
  }
  res_tot <- rbind(res_tot,res)
}
    





library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)
res_varsel1 <- pivot_longer(res_tot, cols = c(1,2,3,4,5,6,7,8,9,10), names_to = "Group", values_to = "Val")
res_varsel1$Group <- factor(res_varsel1$Group, levels = unique(res_varsel1$Group))
res_varsel1$Val = res_varsel1$Val*50 #to report group average

name <- paste("GroupWeights_LinDense","EBBart_WAIC.pdf", sep = "_")
pdf(name, width=7,height=2.5)
bp <- ggplot(res_varsel1, aes(x=Group, y=Val)) +
  geom_boxplot(aes(fill=Model),outlier.shape = NA) + #coord_cartesian(ylim =  c(0.9, 1.3)) +
  #scale_y_continuous(breaks=seq(0.8,2.1,0.2)) +
  theme_light() +
  scale_fill_viridis(discrete = T, option = "B",begin = 0.2, end = 0.7)+
  theme(legend.title = element_blank()) + labs(x="Group",y="Group Weights") +
  geom_hline(yintercept=0.1,linetype=3,linewidth=1)
bp + facet_grid(. ~ Setting)
dev.off()

########### CV vs WAIC ###############

WAIC = results$WAIC
CVPerf = results$CVPerf
Testperf = results$TestPerformance
GroupProbs = results$GroupProbs
remove(results); gc()


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



model <- "flexible"
model <- "rigid"
Setting <- "G = 5, N = 100"



minimaWAIC = apply(WAIC,1, function (x) minimum(x)) # WAIC minimum
minimaCV = apply(CVPerf,1, function (x) minimum(x))

EbcoBART_WAIC <- c()
EBcoBART_CV <- c()
nrep <- length(minimaWAIC)
for (i in 1:nrep){
  EbcoBART_WAIC[i] = Testperf[i,minimaWAIC[i]]
  EBcoBART_CV[i] = Testperf[i,minimaCV[i]]
  #RatioPMSEtest1[i] = PMSEstest[i,minima1[i]]/PMSEstest[i,minima[i]]
}

mean(Testperf[,1]); mean(EbcoBART_WAIC); mean(EBcoBART_CV)

TestRatio = EbcoBART_WAIC/EBcoBART_CV
library(ggplot2); library(viridis)
TestRatio = data.frame(TestRatio = TestRatio)
name <- paste("RatioPerf_WAICvsCV.pdf")
pdf(name, width=4,height=2.5)
bp <- ggplot(TestRatio, aes(y=TestRatio))  +  coord_cartesian(ylim =  c(0.9, 1.1)) +
  geom_boxplot(color ="black", fill = "#69b3a2",fatten =2 ,outlier.shape = NA) +
  scale_x_discrete() +
  #scale_y_continuous(breaks=seq(0.8,2.1,0.2)) +
  theme_light() +
  #scale_fill_viridis(discrete = T, option = "G",begin = 0.4, end = 1)+
  labs(y="PMSE Ratio") 
bp
dev.off()
mean(TestRatio$TestRatio)

getwd()


nrep = 100; ncol = 5

MinGroupProbsCV = matrix(NA, nrow = nrep, ncol = ncol)
for (i in 1:nrep){
  MinGroupProbsCV[i,] <- GroupProbs[[i]][minimaCV[i],]
}

GroupProbs[[1]]


MinGroupProbsWAIC = matrix(NA, nrow = nrep, ncol = ncol)
for (i in 1:nrep){
  MinGroupProbsWAIC[i,] <- GroupProbs[[i]][minimaWAIC[i],]
}

diff = MinGroupProbsWAIC-MinGroupProbsCV

diff = data.frame(diff)
names(diff) = c("1","2","3","4","5")

diff1 <- pivot_longer(diff, cols = c(1,2,3,4,5), names_to = "Group", values_to = "Val")

length(which(abs(diff1$Val)>0.2))
abs(-2)
min(diff[,1])

name <- paste("GroupProbDiff_WAICvsCV.pdf")
pdf(name, width=7,height=2.5)
bp <- ggplot(diff1, aes(x=Group, y=Val)) +
  geom_boxplot(aes(fill=Group),outlier.shape = NA) + coord_cartesian(ylim =  c(-0.25, 0.2)) +
  #scale_y_continuous(breaks=seq(0.8,2.1,0.2)) +
  theme_light() + theme(legend.position = "none") +
  scale_fill_viridis(discrete = T, option = "B",begin = 0.2, end = 0.7)+
  labs(x="Group",y="Diff. Group Weights") 
bp
dev.off()
