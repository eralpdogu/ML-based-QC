
#' A function to train random forest or neural network classifiers for QC data
#'
#' @param data comma-separated (.csv), metric file. It should contain a "peptide" column and the metrics columns. It should also include "Annotations" for each run.
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
#' QcClassifierTrainr(guide.set = sampleData[1:20,], peptide = "VLVLDTDYK", method = "randomforest")

QcClassifierTrain <- function(guide.set, peptide, method){
  
  if(is.null(guide.set))
    return()
  if(!is.data.frame(guide.set)){
    stop(guide.set)
  }
  ###########Simulation#############################################################################
  n<-10000 #incontrol observations 
  Data<-c()
  Data0<-c()
  Data1<-c()
  Data2<-c()
  Data3<-c()
  Data4<-c()
  S0<-c()
  
  #one peptide LVNELTEFAK 
  mean <-c(with(data=guide.set,tapply(RT,INDEX=peptide,FUN=mean))[4],
           with(data=guide.set,tapply(TotalArea,INDEX=peptide,FUN=mean)) [4],
           with(data=guide.set,tapply(MassAccu,INDEX=peptide,FUN=mean)) [4],
           with(data=guide.set,tapply(FWHM,INDEX=peptide,FUN=mean)) [4]
  )
  covar<-cov(guide.set[Guide.set$peptide=="LVNELTEFAK",c(3,6,7,8)])
  #generate in-control observations
  
  S0<-data.frame(idfile=1:(4*n),peptide=rep("LVNELTEFAK",n),mvrnorm(n, mean, covar))
  colnames(S0)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  S0<- reshape(S0, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("GO")
  S0 <- cbind(S0,RESPONSE)
  
  #generate out-of-control observations
  #Logarithmic drift
  Data1<-data.frame(idfile=((1):(n)),peptide=rep("LVNELTEFAK",n),
                    mvrnorm(n, mean+c(1.0*sqrt(covar[1,1]),3.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),1.0*sqrt(covar[4,4])), 
                            covar))
  # Data1<-as.data.frame(matrix(0,n,6))
  # for(i in 1:n){
  #   Data1[i,]<-c(idfile=i,peptide=rep("LVNELTEFAK",1),
  #                mvrnorm(1, mean+c(1.0*sqrt(covar[1,1]),1.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),1.0*sqrt(covar[4,4])), 
  #                        covar))
  # }
  colnames(Data1)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  
  #generate out-of-control observations for a 3 sigma step shift in Total.area---large shift
  Data2<-data.frame(idfile=((n+1):(2*n)),peptide=rep("LVNELTEFAK",n),
                    mvrnorm(n, mean+c(1.0*sqrt(covar[1,1]),3.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),1.0*sqrt(covar[4,4])), 
                            covar))
  colnames(Data2)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  
  #generate out-of-control observations for a 3 sigma fluctuation in Mass.accu---large shift
  covar1<-covar
  covar1[3,3]<-3*covar[3,3]
  Data3<-data.frame(idfile=((2*n+1):(3*n)),peptide=rep("LVNELTEFAK",n),
                    mvrnorm(n, mean, 
                            covar1))
  colnames(Data3)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  
  #generate out-of-control observations for a 1.5 sigma step shift in RT---small shift
  Data4<-data.frame(idfile=(3*n+1):(4*n),peptide=rep("LVNELTEFAK",n),
                    mvrnorm(n, mean+c(1.5*sqrt(covar[1,1]),1.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),1.0*sqrt(covar[4,4])), 
                            covar))
  colnames(Data4)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  
  #Merge all four type of disturbances + in-control observations
  Data0<-rbind(Data1,Data2,Data3,Data4)
  Data0<-reshape(Data0, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("NOGO")
  Data0 <- cbind(Data0,RESPONSE)
  Data0<-rbind(S0,Data0)
  
  #Sample 1000 observations for training
  Train.set<-sample_n(Data0, 10000)
  #Sample 100 observations for test
  Test.set<-sample_n(Data0,1000)
  
  #############Classification########################################################################
  if(method="randomforest"){
  
  #RF model
    fit <- train(RESPONSE ~ ., 
<<<<<<< HEAD
                 data = train_sim_all[,-c(1,2)],  
=======
                 data = Train.set[,-c(1,2)], 
>>>>>>> 166bf69825a99dfe62cd74fb3f586e39a72605ed
                 method="rf",
                 preProcess = c("center", "scale", "nzv"),
                 tuneGrid = data.frame(mtry = 6)) # change mtry as square root of number of predictors
    print(summary(fit))
    
    #Predict RF probabilities and measure performance
    Predict<-predict(fit, Train.set[,-c(1,2,7)]) # change the features here 
    Predict.prob<-predict(fit, Test.set, type="prob")
    print("Test Data",confusionMatrix(Test.set$RESPONSE, Predict))
    
    #Model agnostics
    importance<-varImp(fit)$importance
    plot(varImp(fit))
    
    return(fit) # return the model object 
  }
  
  else if(method="neuralnetwork"){
    #ANN Model
    #Building Neural Network
    library(h2o)
    #generate same set of random numbers (for reproducibility)
    set.seed(121)
    
    #launch h2o cluster
    localH2O <- h2o.init(nthreads = -1)
    
    
    #import r objects to h2o cloud
    train_h2o <- as.h2o(Train.set)
    test_h2o <- as.h2o(Test.set)
    
    set.seed(100)
    
    dl_model_sim1 <- h2o.deeplearning(
      model_id="dl_model_first", 
      training_frame=train_h2o, 
      validation_frame = test_h2o,
      x= colnames(train_h2o[,3:6]),
      y= "RESPONSE",
      activatio="Tanh",  
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


