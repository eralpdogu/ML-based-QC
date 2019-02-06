
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

QCClassifierTrain<- function(guide.set, sim.size){
  
source("auto_add_features.R")
source("sample_density_function.R")
source("boxcox_transformation.R")
source("robust_scaling.R")
source("QcClassifier_data.R")

#function inputs
nmetric<-ncol(guide.set)-2
factor.names = colnames(guide.set[,3:ncol(guide.set)])
#optional 
guide.set$peptide <-as.factor(guide.set$peptide)

#optional 
peptide.colname <-"peptide"

d1 <- QcClassifier_data_step(guide.set,nmetric,factor.names,sim.size*1,peptide.colname)
d2 <- QcClassifier_data_var(guide.set,nmetric,factor.names,sim.size*1,peptide.colname)
d3 <- QcClassifier_data_linear(guide.set,nmetric,factor.names,sim.size*1,peptide.colname)

d<-rbind(d1,d2, d3)
## 80% of the sample size
smp_size <- floor(0.8 * nrow(d))

set.seed(123)
train_ind <- sample(seq_len(nrow(d)), size = smp_size)

train <- d[train_ind,]
test <- d[-train_ind,]


#launch h2o cluster
localH2O <- h2o.init(nthreads = -1)

#import r objects to h2o cloud
train_h2o <- as.h2o(train)
test_h2o <- as.h2o(test)

## run our first predictive model
rf_model <- h2o.randomForest(         ## h2o.randomForest function
  training_frame = train_h2o,        ## the H2O frame for training
  validation_frame = test_h2o,      ## the H2O frame for validation (not required)
  x= colnames(train_h2o),
  y= "RESPONSE",
  model_id = "rf_model",    ## name the model in H2O
  nfolds = 10,
  ntrees = 50,                  ##   not required, but helps use Flow
  ## use a maximum of 200 trees to create the
  ##  random forest model. The default is 50.
  stopping_rounds = 2,           ## Stop fitting new trees when the 2-tree
  ##  average is within 0.001 (default) of 
  ##  the prior two 2-tree averages.
  ##  Can be thought of as a convergence setting
  score_each_iteration = T,
  max_depth = 50,
  seed = 123) 

explainer <- lime(train[,-11], rf_model, n_bins = 5)

results.model<-list()
results.model<-list(train, rf_model, explainer)

return(results.model)
}

#SIMULATION for PERFORMANCE####################################
sequence<-seq(10, 2500, 50)
results<-matrix(NA, 100000, 6)
for(i in sequence){
  cf<- data.frame(h2o.confusionMatrix(QCClassifierTrain(guide.set,i),valid = T),stringsAsFactors = F)
  sens<-cf[1,1]/(cf[1,1]+cf[2,1])
  spec<-cf[2,2]/(cf[1,2]+cf[2,2])
  acc<-(cf[1,1]+cf[2,2])/(cf[1,1]+cf[2,1]+cf[1,2]+cf[2,2])
  err1<-cf[1,3]
  err2<-cf[2,3]
  err3<-cf[3,3]
  results[i,]<-cbind(sens,spec,acc,err1,err2,err3)
}
results<-as.data.frame(cbind(sequence,results[complete.cases(results),]))
results<-results[,c(1,4:6)]
colnames(results)<-c("Simulation.size", "Accuracy", "False positive rate", "False negative rate")

results<-results.DDA
results<-results.SRM

results_melt <- melt(results,id.vars ="Simulation.size")
ggplot(results_melt, aes(Simulation.size, value)) + 
  geom_point()+
  geom_smooth()+
  ylab("Probability")+
  xlab("Simulation size")+
  facet_wrap(~variable,scales = "free")+
  theme_light()

#########################################################################

#SRM example
QCClassifierTrain(guide.set[1:2556,], sim.size=100)
RESPONSE<-NA
Test.set<-guide.set[1:2556,]
Test.set<-cbind(Test.set, RESPONSE)
QCClassifierTest(Test.set)

#DDA example
QCClassifierTrain(guide.set,1000)
QCClassifierTest(rbind(guide.set[101:838,]))
RESPONSE<-NA
test.set<-test.set[,-7]
new.test<-rbind(guide.set[101:838,],test.set[1:100,])
new.test<-cbind(new.test, RESPONSE)
QCClassifierTest(new.test)

g1<-ggplot(test.set[1:500,], aes(idfile, RT)) + 
  geom_line()+
  #ylim(600, 1000)+
  facet_wrap(~peptide,scales = "free")
#QCClassifierInterpret(train, Test.set.scale, rf_model, 48)

ggplot(guide.set, aes(idfile, RT)) + 
  geom_line()+
  ylab("Retention time (sec)")+
  xlab("Run ID")+
  #ylim(600, 1000)+
  facet_wrap(~peptide,scales = "free")+
  colour=variable

ggplot(guide.set, aes(idfile, TotalArea)) + 
  geom_line()+
  ylab("Total peak area")+
  xlab("Run ID")+
  #ylim(600, 1000)+
  facet_wrap(~peptide,scales = "free")+
  theme_light()

ggplot(guide.set, aes(idfile, FWHM)) + 
  geom_line()+
  #ylim(600, 1000)+
  facet_wrap(~peptide,scales = "free")+
  

ggplot(guide.set, aes(idfile, MassAccu)) + 
  geom_line()+
  #ylim(600, 1000)+
  facet_wrap(~peptide,scales = "free")
grid.arrange(g1,g2, ncol=1)


