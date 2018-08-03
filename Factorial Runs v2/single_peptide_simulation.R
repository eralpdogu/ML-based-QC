library(readxl)
library(h2o)
library(caret)
library(MASS)

factorial <- read_xlsx("Factorialcombinatins.xlsx",sheet = 1)


source("auto_add_features.R")

sim.size = 100
tag_neg <- 0
data <- data.frame(NULL)
for(i in 2:nrow(factorial)){
  k = 4 # LVN 
  if(i == 2){
    ###### In cntrol observation ~ 5* sim size  the of the actual 
      sample_data_k <- sample_density(guide.set.scale,levels(guide.set.scale$peptide)[k], sim.size*15)
      colnames(sample_data_k)<- paste(levels(guide.set.scale$peptide)[k],colnames(sample_data_k),sep = ".")
      sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                             sample_data_k,
                             RESPONSE= c("PASS"))
    
  }
  else{
    sample_data<- list()
    ###### In cntrol observation ~ 5* sim size  the of the actual 
      sample_data_k <- sample_density(guide.set.scale,levels(guide.set.scale$peptide)[k], sim.size)
      colnames(sample_data_k)<- paste(levels(guide.set.scale$peptide)[k],colnames(sample_data_k),sep = ".")
      sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                             sample_data_k,
                             RESPONSE= c("PASS"))
    
    
    for(j in 2:5){
      #change in RT Drift for some peptides
      if(factorial[i,j]== "+" & colnames(factorial[i,j])=="RT drift"){ 
            sample_data_k[[paste(sample_data_k$peptide[1],"RT",sep = ".")]] <- 
              sample_data_k[[paste(sample_data_k$peptide[1],"RT",sep = ".")]] +
              runif(1,-2,2)*IQR(sample_data_k[[paste(sample_data_k$peptide[1],"RT",sep = ".")]])
            tag_neg <- 1 
          
      }
      
      #change in Total Area Drift for some peptides
      if(factorial[i,j]== "+" & colnames(factorial[i,j])=="Total Area drift"){ 
            sample_data_k[[paste(sample_data_k$peptide[1],"TotalArea",sep = ".")]] <- 
              sample_data_k[[paste(sample_data_k$peptide[1],"TotalArea",sep = ".")]] +
              runif(1,-2,2)*IQR(sample_data_k[[paste(sample_data_k$peptide[1],"TotalArea",sep = ".")]])
            tag_neg <- 1 
          
          
        
      }
      
      #change in Mass Accu Drift for some peptides
      if(factorial[i,j]== "+" & colnames(factorial[i,j])=="Mass Accu drift"){ 
            sample_data_k[[paste(sample_data_k$peptide[1],"MassAccu",sep = ".")]] <- 
              sample_data_k[[paste(sample_data_k$peptide[1],"MassAccu",sep = ".")]] +
              runif(1,-2,2)*IQR(sample_data_k[[paste(sample_data_k$peptide[1],"MassAccu",sep = ".")]])
            tag_neg <- 1 
      }
      
      #change in FWHM Drift for some peptides
      if(factorial[i,j]== "+" & colnames(factorial[i,j])=="FWHM drift"){ 
            sample_data_k[[paste(sample_data_k$peptide[1],"FWHM",sep = ".")]] <- 
              sample_data_k[[paste(sample_data_k$peptide[1],"FWHM",sep = ".")]] +
              runif(1,-2,2)*IQR(sample_data_k[[paste(sample_data_k$peptide[1],"FWHM",sep = ".")]])
            tag_neg <- 1 
      }
    }# column ends 
  }
  
  data.set <- add_features(sample_data_k)
  #data.set <- data.set[, !duplicated(colnames(data.set))]
  data.set[,"peptide"] <- NULL 
  if(tag_neg == 1){
    data.set[,"RESPONSE"] = rep("FAIL",sim.size)
  }
  tag_neg <- 0
  data <-rbind(data,data.set)
  
}

data <- data[sample(nrow(data), nrow(data)), ] # shuffle the data


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
  x= colnames(train_h2o[,c(1:4,5:21)]),
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

cf<- data.frame(h2o.confusionMatrix(rf_model,valid = T),stringsAsFactors = F)


