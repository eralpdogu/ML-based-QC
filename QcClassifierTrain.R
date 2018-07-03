
#' A function to train random forest or neural network classifiers for QC data
#'
#' @param data comma-separated (.csv), metric file. It should contain a "peptide" column and the metrics columns. It should also include "Annotations" for each run.
#' @param peptide the name of peptide of interest.
#' @param method the method used to model. Two values can be assigned, "randomforest" or "neuralnetwork".
#' @export
#' @import caret pdp ggplot2 MASS dplyr lime h2o 
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- MSstatsQC::DataProcess(S9Site54)
#' head(sampleData)
#' # Find the name of the peptides
#' levels(sampleData$Precursor)
#' # Calculate change point statistics
#' #QcClassifierTrain(guide.set = sampleData,"LVNELTEFAK",method = "neuralnetwork",all_features = F)
QcClassifierTrain <- function(guide.set, peptide,method,all_features, sim.size){
  
  if(is.null(guide.set))
    return()
  if(!is.data.frame(guide.set)){
    stop(guide.set)
  }
  
  guide.set$peptide <-as.factor(guide.set$peptide)
  guide.set.scale<-cbind(guide.set[,1:2],scale(guide.set[,-c(1:2)]))
  j<-which(levels(guide.set.scale$peptide)==peptide)
  
  ###########Simulation#############################################################################
  Data.set<-simulate_disturbances(guide.set.scale, sim.size=100)
  
  ###Splitting Test & Train Data #######################################################################
  ## 75% of the sample size
  smp_size <- floor(0.8 * nrow(Data.set))
  
  ## set the seed to make your partition reproducible
  set.seed(123)
  train_ind <- sample(seq_len(nrow(Data.set)), size = smp_size)
  
  train <- Data.set[train_ind,]
  train<-train[ , order(names(train))]
  
  test <- Data.set[-train_ind,]
  test[ , order(names(test))]
  
  #############Classification########################################################################
  if(method=="randomforest" & all_features == T){
  #RF model
    print("Random Forest: Train Data using all features :")
    fit_all <- train(y=train[,"RESPONSE"],x=subset(train,select = -c(RESPONSE,idfile)), 
                 method="rf",
                 preProcess = c("center", "scale", "nzv"),
                 tuneGrid = data.frame( mtry=floor(sqrt(ncol(train))) )) # change mtry as square root of number of predictors
  
    print(fit_all$results)
  #Model agnostics
    print("Variable Importance")
    plot(varImp(fit_all))
  
    Predict<-predict(fit_all, test)
    Predict.prob<-predict(fit_all, test, type="prob")
    print("Random Forest: Test Data using all features :")
    confusionMatrix(as.factor(test$RESPONSE), Predict,positive='FAIL')
    
    explainer <- lime(subset(train,select = -c(RESPONSE, idfile)), fit_all)
    explanation <- explain(subset(test,select = -c(RESPONSE, idfile)), explainer, n_labels = 1, n_features = 2)
    plot_features(explanation[75:80,])
    
  }
  else if(method=="randomforest" & all_features == F){
    #RF model
    print("Random Forest: Train Data using limited features :")
    fit <- train(y=train[,"RESPONSE"],x=train[2:5], 
                     method="rf",
                     preProcess = c("center", "scale", "nzv"),
                     tuneGrid = data.frame( mtry=floor(sqrt(ncol(train))) )) # change mtry as square root of number of predictors
    
    print(fit$results)
    
    #Model agnostics
    print("Variable Importance")
    plot(varImp(fit))
    
    Predict<-predict(fit, test)
    Predict.prob<-predict(fit, test, type="prob")
    print("Random Forest: Test Data using limited features :")
    confusionMatrix(as.factor(test$RESPONSE), Predict,positive='FAIL')
    
    explainer <- lime(Train.set[2:5], fit)
    explanation <- explain(Test.set[1:4], explainer, n_labels = 1, n_features = 2)
    plot_features(explanation[75:80,])
    
  }
  

  else if(method=="neuralnetwork" & all_features == T){
    #ANN Model
    #generate same set of random numbers (for reproducibility)
    set.seed(121)
    
    #launch h2o cluster
    localH2O <- h2o.init(nthreads = -1)
    
    #import r objects to h2o cloud
    train_h2o <- as.h2o(train)
    test_h2o <- as.h2o(test)
    
    dl_model_sim1 <- h2o.deeplearning(
      model_id="dl_model_first", 
      training_frame=train_h2o, 
      validation_frame = test_h2o,
      x= colnames(train_h2o[,c(2:5,7:9)]),
      y= "RESPONSE",
      activation ="Tanh",  
      hidden=c(20,20), 
      standardize = TRUE, #standardizes the data
      loss= "CrossEntropy",
      stopping_metric="logloss",
      stopping_rounds = 10,
      stopping_tolerance=0.00001,
      adaptive_rate = TRUE,
      shuffle_training_data = TRUE, 
      rate = 0.005, # Defaults to 0.005 adaptive enabled so cannot specify the learning rare 
      mini_batch_size = 1,# defaults to 1 
      epochs=200,
      seed = 123, # give seed 
      export_weights_and_biases = T, # export weights and biases defaults to false
      reproducible = T # Force reproducibility on small data (will be slow - only uses 1 thread). Defaults to FALSE.
    )
    summary(dl_model_sim1)
    plot(dl_model_sim1)
    h2o.shutdown()
    
    return(dl_model_sim1) # return the model object 
    
  }
  else if(method=="neuralnetwork" & all_features == F){
    #ANN Model
    #generate same set of random numbers (for reproducibility)
    set.seed(121)
    
    #launch h2o cluster
    localH2O <- h2o.init(nthreads = -1)
    
    
    #import r objects to h2o cloud
    train_h2o <- as.h2o(train)
    test_h2o <- as.h2o(test)
    
    dl_model_sim1 <- h2o.deeplearning(
      model_id="dl_model_first", 
      training_frame=train_h2o, 
      validation_frame = test_h2o,
      x= colnames(train_h2o[,2:5]),
      y= "RESPONSE",
      activation ="Tanh",  
      hidden=c(20,20), 
      standardize = TRUE, #standardizes the data
      loss= "CrossEntropy",
      stopping_metric="logloss",
      stopping_rounds = 10,
      stopping_tolerance=0.00001,
      adaptive_rate = TRUE,
      shuffle_training_data = TRUE, 
      rate = 0.005, # Defaults to 0.005 adaptive enabled so cannot specify the learning rare 
      mini_batch_size = 1,# defaults to 1 
      epochs=200,
      seed = 123, # give seed 
      export_weights_and_biases = T, # export weights and biases defaults to false
      reproducible = T # Force reproducibility on small data (will be slow - only uses 1 thread). Defaults to FALSE.
    )
    summary(dl_model_sim1)
    plot(dl_model_sim1)
    h2o.shutdown()
    
    return(dl_model_sim1) # return the model object 
    
  }
  
  else{
      
      Print("Illegal Response")
  }

}




