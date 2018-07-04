
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
  setwd("/Users/ed/Dropbox/2. MSstatsQC Paper 2-QCloud")
  Test.set <- read.csv('test_lumos_QCloud_DDA_paper.csv')
  colnames(Test.set)[2]<-"peptide"
  #Test.set<-Test.set[complete.cases(Test.set),]
  
  Test.set$peptide<-as.factor(Test.set$peptide)
  Test.set.scale<-cbind(Test.set[,c(1,2)],scale(Test.set[,-c(1:2)]))
  
  temp.Data<-list()
  for (j in 1:nlevels(Test.set.scale$peptide)){
    temp.Data[[j]]<-Test.set.scale[Test.set.scale$peptide==levels(Test.set.scale$peptide)[j],]
    temp.Data[[j]]<- reshape(temp.Data[[j]], idvar = "idfile", timevar = "peptide", direction = "wide")
  }
  
  Data.set<-add_features(Test.set.scale,temp.Data)
  Test.set<-cbind(Data.set[[1]], 
                  subset(Data.set[[2]],select = -c(idfile)), 
                  subset(Data.set[[3]],select = -c(idfile)),
                  subset(Data.set[[4]],select = -c(idfile)),
                  subset(Data.set[[5]],select = -c(idfile)),
                  subset(Data.set[[6]],select = -c(idfile)),
                  subset(Data.set[[7]],select = -c(idfile)),
                  subset(Data.set[[8]],select = -c(idfile)))
  
  Predict<-predict(fit_all, Test.set, type="prob")
  Predict<-predict(fit_all, Test.set)
  
  explainer <- lime(train, fit_all)
  explanation <- explain(subset(Test.set,select = -c(idfile)), explainer, n_labels = 1, n_features = 2)
  plot_features(explanation)
  

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


