
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
  Results<-list()
  Results_annotated<-list()
  
  for(i in 1:nlevels(Test.set$peptide)){

  Test.set.scale <- Test.set[Test.set$peptide==levels(Test.set$peptide)[i],c(3:(ncol(Test.set)-1))]
  
  guide.set.new<-guide.set[guide.set$peptide==levels(guide.set$peptide)[i],c(3:(ncol(guide.set)))]
  
  for(k in 1:ncol(Test.set.scale)){
  Test.set.scale[,k]=(Test.set.scale[,k]-median(guide.set.new[,k]))/mad(guide.set.new[,k])
  }
  
  #Test.set.scale <- robust.scale(Test.set.scale)
  
  for(k in 1:ncol(Test.set.scale)){Test.set.scale[,k] <- bctrans(Test.set.scale[,k])}
  
  names(Test.set.scale) <- colnames(Test.set[,3:(ncol(Test.set)-1)])

  Test.set.scale <- add_features(Test.set.scale)
  
  Test.set.scale <- Test.set.scale[,order(names(Test.set.scale))]
  
  Test.set.scale.h2o <- as.h2o(Test.set.scale)
  Predict<-as.data.frame(h2o.predict(rf_model, Test.set.scale.h2o, type="prob"))
  Results[[i]]<-Predict$FAIL
  Results_annotated[[i]]<-Predict$predict
  #colnames(Results)[i]<-levels(Test.set$peptide)[i]
  #colnames(Results_annotated)[i]<-levels(Test.set$peptide)[i]
  
  }
  
  Results<-t(plyr::ldply(Results, rbind))
  Results_annotated<-t(plyr::ldply(Results_annotated, rbind))
  colnames(Results)<-levels(Test.set$peptide)
  colnames(Results_annotated)<-levels(Test.set$peptide)

  Test.set.feautures<-cbind(Test.set.scale, Time=1:(dim(Test.set.scale)[1]))
  Test.set.feautures<-melt(Test.set.feautures, id.vars = "Time")
  colnames(Test.set.feautures)[2]<-"Attributes"
   
  g0<-ggplot(Test.set.feautures[-1,], aes(Time, Attributes)) + 
    geom_tile(aes(fill = value), colour = "white") +
    labs(x = "Time",y = "Features and metrics")+
    removeGrid()+
    scale_y_discrete(expand=c(0,0))+
    scale_fill_gradient(low = "white", high = "red",name = "Values")+
    ggtitle(label = "Root causes map")+
    theme(legend.position="bottom")
  
  Results<-data.frame(RUN=1:(dim(Results)[1]), Results)
  # Results_annotated<-data.frame(RUN=1:(dim(Results)[1]), Results_annotated)
  # Results_melt <- melt(Results_annotated,id.vars ="RUN")
  # colors <- c("red","blue")
  # g1<-ggplot(Results_melt, aes(RUN, variable)) + 
  #   geom_tile(aes(fill = value), colour = "white") +
  #   labs(x = "Time",y = "Probability of fail")+
  #   coord_equal()+
  #   ylab("Overall")+
  #   rotateTextX()+
  #   scale_fill_manual(values=colors, name="Label")
  
  Results_melt <- melt(Results[-1,],id.vars ="RUN")
  g2<-ggplot(Results_melt, aes(RUN, variable)) + 
    geom_tile(aes(fill = value), colour = "white") +
    labs(x = "Time",y = "Probability of fail")+
    removeGrid()+
    scale_y_discrete(expand=c(0,0))+
    scale_fill_gradient(low = "white", high = "red",name = "Probability")+
    ggtitle(label = "Decision map")+
    theme(legend.position="bottom")
  gA <- ggplotGrob(g2)
  gB <- ggplotGrob(g0)
  maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
  gA$widths[2:5] <- as.list(maxWidth)
  gB$widths[2:5] <- as.list(maxWidth)
  grid.arrange(gA, ncol=1)  
  # results.test<-list()
  # 
  # explanation_caret <- explain(
  #   x = Test.set.scale, 
  #   explainer = results.model[[3]], 
  #   n_labels = 1,
  #   n_features = 5
  # )
  # results.test<-list(Results, Results_annotated, explanation_caret)
  # return(results.test)
}

  
