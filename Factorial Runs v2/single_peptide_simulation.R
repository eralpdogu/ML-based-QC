
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
new_data <- guide.set
nmetric<-ncol(new_data)-2
factor.names = colnames(new_data[,3:ncol(new_data)])
sim.size = sim.size
#optional 
new_data$peptide <-as.factor(new_data$peptide)

#optional 
peptide.colname <-"peptide"

d1 <- QcClassifier_data_step(new_data,nmetric,factor.names,sim.size,peptide.colname)

#d2 <- QcClassifier_data_var(new_data,nmetric,factor.names,sim.size,peptide.colname)
#d3 <- QcClassifier_data_linear(new_data,nmetric,factor.names,sim.size,peptide.colname)

d<-d1
#d<-rbind(d1,d2,d3)
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
  ntrees = 200,                  ##   not required, but helps use Flow
  ## use a maximum of 200 trees to create the
  ##  random forest model. The default is 50.
  stopping_rounds = 2,           ## Stop fitting new trees when the 2-tree
  ##  average is within 0.001 (default) of 
  ##  the prior two 2-tree averages.
  ##  Can be thought of as a convergence setting
  score_each_iteration = T,     
  seed = 123) 

summary(rf_model)    

boxplot(train,horizontal = T, las=1, cex.axis = 0.5)

cf<- data.frame(h2o.confusionMatrix(rf_model,valid = T),stringsAsFactors = F)
cf

}

QCClassifierTrain(guide.set, 500)
QCClassifierTest(Test.set)
QCClassifierInterpret(train, Test.set.scale, rf_model, 48)
