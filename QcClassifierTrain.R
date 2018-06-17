
#' A function to train random forest or neural network classifiers for QC data
#'
#' @param data comma-separated (.csv), metric file. It should contain a "peptide" column and the metrics columns. It should also include "Annotations" for each run.
#' @param peptide the name of peptide of interest.
#' @param method the method used to model. Two values can be assigned, "randomforest" or "neuralnetwork".
#' @export
#' @import caret pdp ggplot2 MASS dplyr lime
#' @import h2o 
#' @examples
#' # First process the data to make sure it's ready to use
#' sampleData <- MSstatsQC::DataProcess(S9Site54)
#' head(sampleData)
#' # Find the name of the peptides
#' levels(sampleData$Precursor)
#' # Calculate change point statistics
#' #QcClassifierTrain(guide.set = sampleData,"LVNELTEFAK",method = "neuralnetwork",all_features = F)
QcClassifierTrain <- function(guide.set, peptide,method,all_features){
  
  if(is.null(guide.set))
    return()
  if(!is.data.frame(guide.set)){
    stop(guide.set)
  }
  
  guide.set$Precursor <-as.factor(guide.set$Precursor)
  names(guide.set)[6:9] <-c("RT","TotalArea","FWHM","MassAccu")
  guide.set.scale<-cbind(guide.set[,1:2],scale(guide.set[,-c(1:2)]))
  
  
  ###########Simulation#############################################################################
  n<-1000#incontrol observations 
  Data.set<-c()
  Data0<-c()
  Data1<-c()
  Data2<-c()
  Data3<-c()
  Data4<-c()
  S0<-c()
  
  #Peptide :- Name of the peptide
  #convert input peptide to character (input without quotes)
  #peptide = as.character(substitute(Peptide))
  #print(peptide)
  
  #generate in-control observations
  source("sample_density_function.R")
  sample_data <- sample_density(guide.set.scale,peptide, n)
  
  
  
  
  #S0<-data.frame(idfile=1:(4*n),peptide=rep("LVNELTEFAK",n),mvrnorm(n, mean, covar))
  S0<-data.frame(idfile=1:(4*n),peptide=rep(peptide,n),
                 sample_data$sim.sample.RT, sample_data$sim.sample.TotalArea, sample_data$sim.sample.MassAccu, sample_data$sim.sample.FWHM)
  colnames(S0)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  S0<- reshape(S0, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("GO")
  S0 <- cbind(S0,RESPONSE)
  
 
  #generate out-of-control observations
  #Logarithmic drift in FWHM
  sample_data <- sample_density(guide.set.scale,peptide, n)
  
  for(i in 1:n){
    Data1<-rbind(Data1,c(i,rep(peptide,1),
                 sample_data$sim.sample.RT[i], sample_data$sim.sample.TotalArea[i], sample_data$sim.sample.MassAccu[i],
                 sample_data$sim.sample.FWHM[i]+log(i,base=2)*IQR(sample_data$sim.sample.FWHM)))
  }
  colnames(Data1)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  Data1<- as.data.frame(Data1,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data1))){ Data1[,i]<-as.numeric(Data1[,i])}
  
  
  #generate out-of-control observations for a 3 sigma step shift in Total.area---large shift
  sample_data <- sample_density(guide.set.scale,peptide, n)
  
  Data2<-data.frame(idfile = c((n+1):(2*n)),peptide = rep(peptide,n),
                    sample_data$sim.sample.RT, sample_data$sim.sample.TotalArea+2.0*IQR(sample_data$sim.sample.TotalArea), 
                    sample_data$sim.sample.MassAccu, sample_data$sim.sample.FWHM,stringsAsFactors = F)

  colnames(Data2)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  
  #generate out-of-control observations for a 3 sigma fluctuation in Mass.accu---large variation, change in shape of the density
  sample_data <- sample_density(guide.set.scale,peptide, n)
  
  Data3<-data.frame(idfile=((2*n+1):(3*n)),peptide=rep(peptide,n),
                    sample_data$sim.sample.RT, sample_data$sim.sample.TotalArea, sample_data$sim.sample.MassAccu^(3), 
                    sample_data$sim.sample.FWHM,stringsAsFactors = F)
  
  colnames(Data3)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  
  #generate out-of-control observations for a 1.5 sigma step shift in RT---small shift
  sample_data <- sample_density(guide.set.scale,peptide, n)
  
  Data4<-data.frame(idfile=(3*n+1):(4*n),peptide=rep(peptide,n),
                    sample_data$sim.sample.RT+0.75*IQR(sample_data$sim.sample.RT), sample_data$sim.sample.TotalArea, 
                    sample_data$sim.sample.MassAccu, sample_data$sim.sample.FWHM,stringsAsFactors = F)
  
  colnames(Data4)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  
  #Merge all four type of disturbances + in-control observations
  Data.set<-rbind(Data1,Data2,Data3,Data4)
  Data.set<-reshape(Data.set, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("NOGO",n),rep("NOGO",n),rep("NOGO",n),rep("NOGO",n))
  Data.set <- cbind(Data.set,RESPONSE)
  Data.set<-rbind(S0,Data.set)
  Data.set$idfile<-1:(8*n)
  
  
  #############Include MR and CUSUMs#################################################################
  #for(j in 1:(nlevels(guide.set$peptide)*nmetrics)){
  for(j in 2:5){
  v <- numeric(length(Data.set[,j]))
  MR <- numeric(length(Data.set[,j]))
  CUSUMpoz.m <- numeric(length(Data.set[,j]))
  CUSUMneg.m <- numeric(length(Data.set[,j]))
  CUSUMpoz.v <- numeric(length(Data.set[,j]))
  CUSUMpoz.m[1]<-0
  CUSUMneg.m[1]<-0
  CUSUMpoz.v[1]<-0
  MR[1]<-0
  v[1]<-0
  k<-0.5
  for(i in 2:length(Data.set[,j])) {
    MR[i] <- abs(Data.set[i,j]-Data.set[(i-1),j])
    v[i] <- (sqrt(abs(Data.set[i,j]))-0.822)/0.349
    CUSUMpoz.m[i] <- max(0,(Data.set[i,j]-(k)+CUSUMpoz.m[i-1]))
    CUSUMneg.m[i] <- max(0,((-k)-Data.set[i,j]+CUSUMneg.m[i-1]))
    CUSUMpoz.v[i] <- max(0,(v[i]-(k)+CUSUMpoz.v[i-1]))
  }
 }
  addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMpoz.v)
  colnames(addfeatures) <- paste(colnames(Data.set[j]),colnames(addfeatures), sep = ".")
  Data.set<-cbind(Data.set,addfeatures) # full dataset with external features
  
  ###Splitting Test & Train Data #######################################################################
  ## 75% of the sample size
  smp_size <- floor(0.8 * nrow(Data.set))
  
  ## set the seed to make your partition reproducible
  set.seed(123)
  train_ind <- sample(seq_len(nrow(Data.set)), size = smp_size)
  
  
  train <- Data.set[train_ind,]
  test <- Data.set[-train_ind,]

  
  #############Classification########################################################################
  if(method=="randomforest" & all_features == T){
  #RF model
    print("Random Forest: Train Data using all features :")
    fit_all <- train(y=train[,"RESPONSE"],x=subset(train,select = -c(RESPONSE)), 
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
    confusionMatrix(as.factor(test$RESPONSE), Predict,positive='NOGO')
    
    explainer <- lime(subset(train,select = -c(RESPONSE)), fit_all)
    explanation <- explain(subset(test,select = -c(RESPONSE)), explainer, n_labels = 1, n_features = 2)
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
    confusionMatrix(as.factor(test$RESPONSE), Predict,positive='NOGO')
    
    explainer <- lime(train[2:5], fit)
    explanation <- explain(test[2:5], explainer, n_labels = 1, n_features = 2)
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




