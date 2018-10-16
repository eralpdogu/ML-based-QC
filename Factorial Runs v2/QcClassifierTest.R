
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

# QcClassifierTest<- function(guide.set, Test.set, peptide, method, nmetrics){
  
  Test.set<-Data.set
  Test.set$peptide<-as.factor(Test.set$peptide)
  
  Results<-as.data.frame(matrix(0,dim(Test.set)[1]/nlevels(Test.set$peptide),nlevels(Test.set$peptide)))
  Results_annotated<-as.data.frame(matrix(0,dim(Test.set)[1]/nlevels(Test.set$peptide),nlevels(Test.set$peptide)))
  
  for(i in 1:nlevels(Test.set$peptide)){
  sample_data_k <- guide.set[guide.set$peptide==levels(guide.set$peptide)[i],3:6]

  Test.set.scale<-Test.set[Test.set$peptide==levels(Test.set$peptide)[i],3:6]
  
  for(j in 1:ncol(Test.set.scale)){
    Test.set.scale[,j]=(Test.set.scale[,j]-median(sample_data_k[,j]))/mad(sample_data_k[,j])
  }
  
  sample_data_k<-robust.scale(sample_data_k)
  
  lambda.fm1 <- boxCox(sample_data_k$RT ~ 1, family ="yjPower", plotit = FALSE)
  lambda.max <- lambda.fm1$x[which.max(lambda.fm1$y)]
  Test.set.scale$RT = yjPower(Test.set.scale$RT, lambda=lambda.max, jacobian.adjusted=FALSE)
  
  lambda.fm1 <- boxCox(sample_data_k$TotalArea ~ 1, family ="yjPower", plotit = FALSE)
  lambda.max <- lambda.fm1$x[which.max(lambda.fm1$y)]
  Test.set.scale$TotalArea = yjPower(Test.set.scale$TotalArea, lambda=lambda.max, jacobian.adjusted=FALSE)
  
  lambda.fm1 <- boxCox(sample_data_k$MassAccu~ 1, family ="yjPower", plotit = FALSE)
  lambda.max <- lambda.fm1$x[which.max(lambda.fm1$y)]
  Test.set.scale$MassAccu = yjPower(Test.set.scale$MassAccu, lambda=lambda.max, jacobian.adjusted=FALSE)
  
  lambda.fm1 <- boxCox(sample_data_k$FWHM ~ 1, family ="yjPower", plotit = FALSE)
  lambda.max <- lambda.fm1$x[which.max(lambda.fm1$y)]
  Test.set.scale$FWHM = yjPower(Test.set.scale$FWHM, lambda=lambda.max, jacobian.adjusted=FALSE)
  
  names(Test.set.scale) <- c("RT", "TotalArea", "MassAccu", "FWHM")
  
  Test.set.scale <- add_features(Test.set.scale)
  Test.set.scale.h2o <- as.h2o(Test.set.scale)
  Predict<-as.data.frame(h2o.predict(rf_model, Test.set.scale.h2o, type="prob"))
  Results[,i]<-Predict$FAIL
  Results_annotated[,i]<-Predict$predict
  colnames(Results)[i]<-levels(Test.set$peptide)[i]
  colnames(Results_annotated)[i]<-levels(Test.set$peptide)[i]
  
  #explainer <- lime(train, rf_model)
  #explanation <- explain(Test.set.scale, explainer, n_labels = 1, n_features = 2)
  #plot_features(explanation[48:50,])
  
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

#}


