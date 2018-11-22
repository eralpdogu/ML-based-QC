
library(h2o)
library(MASS)
library(ggExtra)
library(ggplot2)
library(gridExtra)
library(stats)
library(FrF2) 
library(car)

source("auto_add_features.R")
source("sample_density_function.R")
nmetric<-ncol(guide.set)-2

factorial <- FrF2(2^nmetric, nmetric,factor.names=colnames(guide.set[,3:ncol(guide.set)]))

guide.set$peptide <-as.factor(guide.set$peptide)

sim.size = 100
tag_neg <- 0

data <- data.frame(NULL)

for(i in 1:nrow(factorial)){
  data.set <- data.frame(NULL)
  if(factorial[i,]== rep(-1,nmetric)){
    ####### In cntrol observation ~ 5* sim size  the of the actual 
    sample_data_k <- sample_density(guide.set, sim.size)
  }
  
  else{
    ###### Base Data set to begin with 
    sample_data_k <- sample_density(guide.set, sim.size)
    #sample_data_k <- robust.scale(sample_data_k)
    
    for(j in 1:ncol(sample_data_k)){
      #change in a metric for some peptides
      if(factorial[i,j]== "1" & colnames(factorial[i,j])==colnames(sample_data_k)[j]){ 
        beta=runif(sim.size,-5,5)
        sample_data_k[,j] <- sample_data_k[,j] + beta*mad(sample_data_k[,j])
        tag_neg <- 1 
    
    }# column ends 
    }
  }
  data.set <- rbind(add_features(sample_data_k))
  #data.set[,"peptide"] <- NULL 
  if(tag_neg == 1){
    data.set$RESPONSE <- c("FAIL")
    tag_neg <- 0
  }
  else{
    data.set$RESPONSE <- c("PASS")
  }
  data <- data[,order(names(data))]
  data.set <- data.set[,order(names(data.set))]
  data <-rbind(data,data.set)
}

data <- data[sample(nrow(data), nrow(data)), ] # shuffle the data
data$RESPONSE <- as.factor(data$RESPONSE)

## 80% of the sample size
smp_size <- floor(0.8 * nrow(data))

set.seed(123)
train_ind <- sample(seq_len(nrow(data)), size = smp_size)

train <- data[train_ind,]
test <- data[-train_ind,]


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

