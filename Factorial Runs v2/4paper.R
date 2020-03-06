
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

#FIGURE 4-5
#SRM example
guide.set<-guide.set.SRM[1:2000,]
Test.set<-test.set.SRM

MSstatsQC.ML.trainR(guide.set, sim.size=2000)

Test.set<-rbind(guide.set[1:500,],test.set.SRM[1:200,])
MSstatsQC.ML.testR(Test.set, guide.set)

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

MSstatsQC.ML.trainR(guide.set[,], sim.size=1000, guide.set.annotations = NULL)
MSstatsQC.ML.trainR(guide.set, sim.size=1000, guide.set.annotations = guide.set.DDA.anno)
MSstatsQC.ML.trainR(guide.set, sim.size=100)

new.test<-rbind(guide.set[,], test.set[test.set$idfile>41234&test.set$idfile<41505 ,])
new.test<-rbind(guide.set[201:838,], test.set[1:200,])
MSstatsQC.ML.testR(new.test, guide.set)

guide.set %>%
  group_by(peptide) %>%
  summarise_each(funs(mean))

test.set.DDA2 %>%
  group_by(peptide) %>%
  summarise_each(funs(mean))
