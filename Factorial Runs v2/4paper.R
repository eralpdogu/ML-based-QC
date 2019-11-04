
# A -----------------------------------------------------------------------



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
SimData<-data.frame(Data.set[,1:2], Annotations=NA, Data.set[,3:6])
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
  geom_line()+ 
  # geom_smooth(data=filter(Simdata_melt, 
  #                         Simdata_melt$variable == "RetentionTime"&Simdata_melt$Run>25), 
  #             aes(Run, value), method = "lm") +
  geom_point(data=filter(Simdata_melt, Simdata_melt$variable == "RetentionTime"), 
             aes(Run, value))+
  geom_smooth(data=filter(Simdata_melt, Simdata_melt$variable == "RetentionTime"), 
             aes(Run, value), method="glm")+
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
Test.set<-test.set.SRM

MSstatsQC.ML.trainR(guide.set, sim.size=1000)

RESPONSE<-NA
Test.set<-rbind(guide.set[1:500,],test.set.SRM[1:200,])
test.set<-cbind(Test.set, RESPONSE)
MSstatsQC.ML.testR(Test.set, guide.set)

RESPONSE<-NA
Test.set<-rbind(test.set.SRM[1:500,])
test.set<-cbind(Test.set, RESPONSE)
MSstatsQC.ML.testR(Test.set, guide.set)

MSstatsQC.ML.trainR(Test.set, sim.size=1000)

guide.set %>%
  group_by(peptide) %>%
  summarise_each(funs(mean))

test.set.SRM %>%
  group_by(peptide) %>%
  summarise_each(funs(mean))

ggplot(guide.set, aes(RT))+ 
  geom_density()+ 
  facet_wrap(~peptide)

ggplot(test.set.SRM[1:200,], aes(RT))+ 
  geom_density()+ 
  facet_wrap(~peptide)

ggplot(guide.set, aes(TotalArea))+ 
  geom_density()+ 
  facet_wrap(~peptide)

ggplot(test.set.SRM[1:200,], aes(TotalArea))+ 
  geom_density()+ 
  facet_wrap(~peptide)

#DDA example
guide.set<-guide.set.DDA[,]
test.set<-test.set.DDA2[,]

guide.set<-test.set.DDA[,]
RESPONSE<-NA
guide.set<-cbind(guide, RESPONSE)

MSstatsQC.ML.trainR(guide.set, sim.size=1000, guide.set.annotations = NULL)
MSstatsQC.ML.trainR(guide.set, sim.size=1000, guide.set.annotations = guide.set.DDA.anno)
MSstatsQC.ML.trainR(test.set, sim.size=1000)

new.test<-rbind(guide.set[,], test.set[test.set$idfile>41234&test.set$idfile<41505 ,])
new.test<-rbind(guide.set, test.set[1:200,])
MSstatsQC.ML.testR(new.test, guide.set)

new.test<-rbind(guide.set, test.set)
MSstatsQC.ML.testR(new.test, test.set)

guide.set %>%
  group_by(peptide) %>%
  summarise_each(funs(mean))

test.set.DDA2 %>%
  group_by(peptide) %>%
  summarise_each(funs(mean))
