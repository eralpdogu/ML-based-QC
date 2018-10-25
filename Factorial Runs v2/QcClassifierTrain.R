#' A function to train random forest classifiers for QC data
#'
#' @param guide.set comma-separated (.csv), metric file. It should contain a "peptide" column and the metrics columns. It should also include "Annotations" for each run.
#' @param Test.set comma-separated (.csv), metric file. It should contain a "peptide" column and the metrics columns. It should also include "Annotations" for each run.
#' @param peptide the name of peptide of interest.
#' @param method the method used to model. Two values can be assigned, "randomforest" or "neuralnetwork".
#' @export
#' @import caret pdp ggplot2 MASS dplyr
#' @import h2o
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- MSstatsQC::DataProcess(S9Site54)
#' head(sampleData)
#' # Find the name of the peptides
#' levels(sampleData$Precursor)
#' # Calculate change point statistics
#' QcClassifierTrain(guide.set = sampleData[1:20,])

QcClassifierTrain <- function(guide.set, sim.size=1000){
  
  source("auto_add_features.R")
  source("sample_density_function.R")
  
  factorial <- read_xlsx("Factorialcombinatins.xlsx",sheet = 1)
  
  guide.set$peptide <-as.factor(guide.set$peptide)
  
  tag_neg <- 0
  
  data <- data.frame(NULL)
  
  for(i in 2:nrow(factorial)){
    data.set <- data.frame(NULL)
    if(i == 2){
      ####### In cntrol observation ~ 5* sim size  the of the actual 
      sample_data_k <- sample_density(guide.set, sim.size*15)
    }
    
    else{
      ###### Base Data set to begin with 
      sample_data_k <- sample_density(guide.set.scale, sim.size)

      for(j in 2:5){
        #change in RT Drift for some peptides
        if(factorial[i,j]== "+" & colnames(factorial[i,j])=="RT drift"){ 
          beta=runif(sim.size,-3,3)
          sample_data_k$RT <- sample_data_k$RT + beta*mad(sample_data_k$RT)
          tag_neg <- 1 
        }
        
        #change in Total Area Drift for some peptides
        if(factorial[i,j]== "+" & colnames(factorial[i,j])=="Total Area drift"){ 
          beta=runif(sim.size,-3,3)
          sample_data_k$TotalArea <- sample_data_k$TotalArea + beta*mad(sample_data_k$TotalArea)
          tag_neg <- 1 
        }
        
        #change in Mass Accu Drift for some peptides
        if(factorial[i,j]== "+" & colnames(factorial[i,j])=="Mass Accu drift"){ 
          beta=runif(sim.size,-3,3)
          sample_data_k$MassAccu <- sample_data_k$MassAccu + beta*mad(sample_data_k$MassAccu)
          tag_neg <- 1 
        }
        #change in FWHM Drift for some peptides
        if(factorial[i,j]== "+" & colnames(factorial[i,j])=="FWHM drift"){
          beta=runif(sim.size,-3,3)
          sample_data_k$FWHM <- sample_data_k$FWHM + beta*mad(sample_data_k$FWHM)
          tag_neg <- 1 
        }
        
      }# column ends 
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
    x= colnames(train_h2o[,c(1:4,5:17)]),
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