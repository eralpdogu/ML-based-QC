
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
#' QcClassifierTrain(guide.set = sampleData[1:20,], peptide = "LVNELTEFAK", method = "randomforest")

  QCClassifierTest<- function(Test.set){
  
  source("auto_add_features.R")
  source("robust_scaling.R")
  
  Test.set$peptide<-as.factor(Test.set$peptide)
  
  Results<-as.data.frame(matrix(0,max(table(Test.set$peptide)),nlevels(Test.set$peptide)))
  Results_annotated<-as.data.frame(matrix(0,max(table(Test.set$peptide)),nlevels(Test.set$peptide)))
  
  for(i in 1:nlevels(Test.set$peptide)){

  Test.set.scale <- Test.set[Test.set$peptide==levels(Test.set$peptide)[i],c(3:(ncol(Test.set)-1))]
  
  guide.set.new<-guide.set[guide.set$peptide==levels(guide.set$peptide)[i],c(3:(ncol(guide.set)))]
  
  for(k in 1:ncol(Test.set.scale)){
  Test.set.scale[,k]=(Test.set.scale[,k]-median(guide.set.new[,k]))/mad(guide.set.new[,k])
  }
  
  #Test.set.scale <- robust.scale(Test.set.scale)
  
  for(k in 1:ncol(Test.set.scale)){Test.set.scale[,k] <- bctrans(Test.set.scale[,k])}
  
  names(Test.set.scale) <- c("RT", "TotalArea", "MassAccu", "FWHM")

  Test.set.scale <- add_features(Test.set.scale)
  
  Test.set.scale <- Test.set.scale[,order(names(Test.set.scale))]
  
  Test.set.scale.h2o <- as.h2o(Test.set.scale)
  Predict<-as.data.frame(h2o.predict(rf_model, Test.set.scale.h2o, type="prob"))
  Results[,i]<-Predict$FAIL
  Results_annotated[,i]<-Predict$predict
  colnames(Results)[i]<-levels(Test.set$peptide)[i]
  colnames(Results_annotated)[i]<-levels(Test.set$peptide)[i]
  
  }
  
  boxplot(Test.set.scale,horizontal = T, las=1, cex.axis = 0.5)
  
  Results<-data.frame(RUN=1:(dim(Test.set)[1]/nlevels(Test.set$peptide)), Results)
  Results_annotated<-data.frame(RUN=1:(dim(Test.set)[1]/nlevels(Test.set$peptide)), Results_annotated)
  
  Results_melt <- melt(Results_annotated,id.vars ="RUN")
  colors <- c("red","blue")
  g1<-ggplot(Results_melt, aes(RUN, variable)) + 
    geom_tile(aes(fill = value), colour = "white") +
    labs(x = "RUN",y = "Probability of FAIL")+
    coord_equal()+
    ylab("Overall")+
    rotateTextX()+
    scale_fill_manual(values=colors, name="Label")
  
  Results_melt <- melt(Results,id.vars ="RUN")
  g2<-ggplot(Results_melt, aes(RUN, variable)) + 
    geom_tile(aes(fill = value), colour = "white") +
    labs(x = "RUN",y = "Probability of FAIL")+
    rotateTextX()+
    removeGrid()+
    scale_fill_gradient(low = "blue", high = "red",name = "Probability") 
  
  grid.arrange(g1,g2)

}


