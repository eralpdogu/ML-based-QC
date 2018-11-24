
library(h2o)
library(MASS)
library(ggExtra)
library(ggplot2)
library(gridExtra)
library(stats)
library(FrF2) 
library(car)

guide.set1 <- guide.set
guide.set2 <- guide.set[,1:5]

source("auto_add_features.R")
source("sample_density_function.R")
source("boxcox_transformation.R")
source("robust_scaling.R")

#function inputs
new_data <- guide.set2
nmetric<-ncol(new_data)-2
factor.names = colnames(new_data[,3:ncol(new_data)])
sim.size = 100
#optional 
new_data$peptide <-as.factor(new_data$peptide)

QcClassifie_data <- function(data,nmetrics,factor.names,sim.size,peptide.colname){
  
  if(!is.factor(new_data[,paste(peptide.colname)])){
    new_data$peptide <-as.factor(new_data$peptide.colname)  
  }
  #factorial matrix  
  factorial <- FrF2(2^nmetric, nmetric,factor.names=factor.names)
  
  
  tag_neg <- 0
  data <- data.frame(NULL)
  
  for(i in 1:nrow(factorial)){
    data.set <- data.frame(NULL)
    if(all(factorial[i,]== rep(-1,nmetric))){
      ####### In cntrol observation ~ 5* sim size  the of the actual 
      sample_data_k <- sample_density(new_data, sim.size*(2^(nmetric)-1))
    }
    
    else{
      ###### Base Data set to begin with 
      sample_data_k <- sample_density(new_data, sim.size)
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
    data.set <- rbind(data.set,add_features(sample_data_k))
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
  
  return(data) 
}

# Test Cases : 
new_data <- guide.set2
nmetric<-ncol(new_data)-2
factor.names = colnames(new_data[,3:ncol(new_data)])
sim.size = 100
#optional 
peptide.colname <-"peptide"


d <- QcClassifie_data(new_data,nmetric,factor.names,sim.size = 100,peptide.colname)


  

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

