
library(h2o)
library(MASS)
library(ggExtra)
library(ggplot2)
library(gridExtra)
library(stats)
library(FrF2) 
library(car)
library(reshape2)
library(lime)
library(dplyr)

setwd("/Users/ed/GitHub/ML-based-QC/Factorial Runs v2")

#FIGURE 1
#raw data plots
SimData<-cbind(SimData, SimData[,4])
SimData<-SimData[,-4]
colnames(SimData)<-c("Run", "Precursor", "Annotations", "TotalArea", "MassAccuracy", "FWHM", "RetentionTime")

Sim.set.scale <- SimData[SimData$Precursor==levels(SimData$Precursor)[1],c(4:ncol(SimData))]
for(k in 1:ncol(Sim.set.scale)){
  Sim.set.scale[,k]=(Sim.set.scale[,k]-median(Sim.set.scale[,k]))/mad(Sim.set.scale[,k])
}
Sim.set<-Sim.set.scale

for(i in 2:nlevels(SimData$Precursor)){
  
  Sim.set.scale <- SimData[SimData$Precursor==levels(SimData$Precursor)[i],c(4:ncol(SimData))]
  
  for(k in 1:ncol(Sim.set.scale)){
    Sim.set.scale[,k]=(Sim.set.scale[,k]-median(Sim.set.scale[,k]))/mad(Sim.set.scale[,k])
  }
  Sim.set<-rbind(Sim.set, Sim.set.scale)
}
Sim.set<-cbind(SimData[,1:2], Sim.set)


Simdata_melt <- melt(Sim.set,id.vars =c("Precursor","Run"))

ggplot(Simdata_melt, aes(Run, value, color=variable)) + 
  geom_smooth(method="lm")+ 
  # geom_smooth(data=filter(Simdata_melt, 
  #                         Simdata_melt$variable == "RetentionTime"&Simdata_melt$Run>25), 
  #             aes(Run, value), method = "lm") +
  geom_point(data=filter(Simdata_melt, Simdata_melt$variable == "RetentionTime"), 
             aes(Run, value))+
  
  ylab(NULL)+
  xlab("Time")+
  ylim(-3, 3)+
  facet_wrap(~Precursor,scales = "free", ncol = 4)+
  scale_color_manual(values = c("#F0E442", "#0072B2", "#CC79A7", "#D55E00"))+
  labs(color = "Metric")+
  theme(legend.position="bottom", panel.background = element_blank(),
        plot.background = element_blank(), plot.margin = unit(c(0.1,0,0,0), "cm"),
        axis.ticks.length = unit(0, "pt"))

#FIGURE 4-5
#SRM example
guide.set<-guide.set.SRM
MSstatsQC.ML.trainR(guide.set, sim.size=10)

RESPONSE<-NA
Test.set<-rbind(guide.set[1:300,],test.set.SRM[1:200,])
test.set<-cbind(Test.set, RESPONSE)
MSstatsQC.ML.testR(Test.set, guide.set)

RESPONSE<-NA
Test.set<-cbind(test.set.SRM, RESPONSE)
MSstatsQC.ML.testR(test.set.SRM, guide.set)

#DDA example
guide.set<-guide.set.DDA
test.set<-test.set.DDA2
MSstatsQC.ML.trainR(guide.set, sim.size=1000, guide.set.annotations = guide.set.DDA.anno)
MSstatsQC.ML.trainR(guide.set, sim.size=1000)

new.test<-rbind(guide.set[1:200,], test.set[test.set$idfile>41234&test.set$idfile<41505 ,])
new.test<-rbind(guide.set, test.set[1:500,])
MSstatsQC.ML.testR(new.test, guide.set)

new.test<-rbind(guide.set, test.set)
MSstatsQC.ML.testR(new.test, guide.set)

