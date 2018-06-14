
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
#' QcClassifierTrainr(guide.set = sampleData[1:20,], peptide = "LVNELTEFAK", method = "randomforest")

QcClassifierTrain <- function(guide.set, peptide, method, nmetrics){
  
  if(is.null(guide.set))
    return()
  if(!is.data.frame(guide.set)){
    stop(guide.set)
  }
  
  guide.set$peptide<-as.factor(guide.set$peptide)
  ###########Simulation#############################################################################
  n<-1000 #incontrol observations 
  Train.set<-c()
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
  covar<-cov(guide.set[guide.set$peptide=="LVNELTEFAK",c(3,6,7,8)])
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
  Data1<-as.data.frame(matrix(0,n,6))
  for(i in 1:n){
  Data1[i,]<-c(idfile=i,peptide=rep("LVNELTEFAK",1),
                  mvrnorm(1, mean+c(1.0*sqrt(covar[1,1]),1.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),log(i,base=2)*sqrt(covar[4,4])), 
                          covar))
  }
  for (i in 3:ncol(Data1)){ Data1[,i]<-as.numeric(Data1[,i])}
 
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
  Train.set<-rbind(Data1,Data2,Data3,Data4)
  Train.set<-reshape(Train.set, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("NOGO")
  Train.set <- cbind(Train.set,RESPONSE)
  Train.set<-rbind(S0,Train.set)
  
  #############Include MR and CUSUMs#################################################################
  #for(j in 1:(nlevels(guide.set$peptide)*nmetrics)){
  for(j in 2:5){

  v <- numeric(length(Train.set[,j]))
  t <- numeric(length(Train.set[,j]))
  Cpoz.mean <- numeric(length(Train.set[,j]))
  Cneg.mean <- numeric(length(Train.set[,j]))
  Cpoz.var <- numeric(length(Train.set[,j]))
  Cneg.var <- numeric(length(Train.set[,j]))
  
  k<-0.5
  for(i in 2:length(Train.set[,j])) {
    Cpoz.mean[i] <- max(0,(Train.set[i,j]-(k)+Cpoz[i-1]))
    Cneg.mean[i] <- max(0,((-k)-Train.set[i,j]+Cneg[i-1]))
  }
  
    for(i in 2:length(Train.set[,j])) {
      v[i] <- (sqrt(abs(Train.set[i,j]))-0.822)/0.349
    }
  
    for(i in 2:length(Train.set[,j])) {
      Cpoz.var[i] <- max(0,(v[i]-(k)+Cpoz[i-1]))
      Cneg.var[i] <- max(0,((-k)-v[i]+Cneg[i-1]))
    }
  
  for(i in 2:length(Train.set[,j])) {
    t[i] <- abs(Train.set[i,j]-Train.set[(i-1),j])
  }
  addfeatures<-cbind(Cpoz.mean,Cneg.mean,Cpoz.var,Cneg.var,t)
  colnames(addfeatures) <- paste(levels(guide.set$peptide)[4], colnames(addfeatures), sep = ".")
  }
  
  Train.set<-cbind(Train.set,addfeatures)
  #############Classification########################################################################
  if(method="randomforest"){
  
  #RF model
    fit <- train(RESPONSE ~ ., 
                 data = Train.set[,-c(1,2)], 
                 method="rf",
                 preProcess = c("center", "scale", "nzv"),
                 tuneGrid = data.frame(mtry = 6)) # change mtry as square root of number of predictors
    
    #Model agnostics
    importance<-varImp(fit)$importance
    plot(varImp(fit))
    
    return(fit) # return the model object 
  }
  
  else if(method="neuralnetwork"){
    #ANN Model
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


