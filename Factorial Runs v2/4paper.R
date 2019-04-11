
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
MSstatsQC.ML.trainR(guide.set[1:2556,], sim.size=2000)
RESPONSE<-NA
Test.set<-guide.set[1:2556,]
Test.set<-cbind(Test.set, RESPONSE)
MSstatsQC.ML.testR(Test.set, Test.set)

#DDA example
guide.set<-guide.set.DDA
MSstatsQC.ML.trainR(guide.set,1000)
MSstatsQC.ML.testR(rbind(guide.set[1:838,]), guide.set)
RESPONSE<-NA
new.test<-rbind(guide.set[1:838,], test.set[1:160,])
new.test<-cbind(new.test, RESPONSE)
MSstatsQC.ML.testR(new.test, guide.set)

MSstatsQC.ML.trainR(test.set, 1000)
MSstatsQC.ML.testR(test.set[,], test.set)

