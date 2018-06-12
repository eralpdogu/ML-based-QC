
#' A function to test random forest or neural network classifiers for QC data
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
#' QcClassifierTrainr(guide.set = sampleData[1:20,], peptide = "LVNELTEFAK", method = "randomforest")

QcClassifierTest<- function(guide.set, Test.set, peptide, method, nmetrics){
  
  if(is.null(Test.set))
    return()
  if(!is.data.frame(Test.set)){
    stop(Test.set)
  }
  
  Test.set$peptide<-as.factor(Test.set$peptide)
 
  Test.set<-reshape(Test.set, idvar = "idfile", timevar = "peptide", direction = "wide")
  
  #############Classification########################################################################
  if(method="randomforest"){
  
  #RF model
  fit<-QcClassifierTrain(QcClassifierTrain(guide.set, peptide = "LVNELTEFAK", method = "randomforest"))
  
  #Predict RF probabilities and measure performance
  Predict.prob<-predict(fit, Test.set, type="prob")
  }
  
  #NN model
  else if(method="neuralnetwork"){
  fit<-QcClassifierTrain(QcClassifierTrain(guide.set, peptide = "LVNELTEFAK", method = "randomforest"))
  
  #Predict NN probabilities and measure performance
  #Predict.prob<-predict(fit, Test.set, type="prob")  
  }
    else{
      
      Print("Illegal Response")
    }

}


