
#' A function to test random forest classifiers for QC data
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

  MSstatsQC.ML.testR<- function(Test.set, guide.set,address=""){
  
  source("auto_add_features.R")
  source("robust_scaling.R")
  source("boxcox_transformation.R")
  
    if (address != FALSE) {
      allfiles <- list.files()
      
      num <- 0
      filenaming <- paste0(address,"MSstatsQC.ML.Plots")
      finalfile <- paste0(address,"MSstatsQC.ML.Plots.pdf")
      
      while (is.element(finalfile, allfiles)) {
        num <- num + 1
        finalfile <- paste0(paste(filenaming, num, sep="-"), ".pdf")
      }	
      
      pdf(finalfile, width=20, height=20)
    }
    
  Test.set$peptide<-as.factor(Test.set$peptide)
  guide.set$peptide<-as.factor(guide.set$peptide)
  Results<-list()
  Results_annotated<-list()
  Test.set.features<-list()
  interpret.plots<-list()
  
  for(i in 1:nlevels(Test.set$peptide)){

  Test.set.scale <- Test.set[Test.set$peptide==levels(Test.set$peptide)[i],c(3:(ncol(Test.set)-1))]
  
  guide.set.new<-guide.set[guide.set$peptide==levels(guide.set$peptide)[i],c(3:(ncol(guide.set)))]
  
  for(k in 1:ncol(Test.set.scale)){
  Test.set.scale[,k]=(Test.set.scale[,k]-median(guide.set.new[,k]))/mad(guide.set.new[,k])
  }
  
  guide.set.new <- robust.scale(guide.set.new)
  
  for(k in 1:ncol(Test.set.scale)){Test.set.scale[,k] <- bctrans.test((guide.set.new[,k]),Test.set.scale[,k])}
  
  names(Test.set.scale) <- colnames(Test.set[,3:(ncol(Test.set)-1)])

  Test.set.scale <- add_features(Test.set.scale)
  
  Test.set.scale <- Test.set.scale[,order(names(Test.set.scale))]
  
  Test.set.scale.h2o <- as.h2o(Test.set.scale)
  Predict<-as.data.frame(h2o.predict(results.model[[2]], Test.set.scale.h2o, type="prob"))
  Results[[i]]<-Predict$FAIL
  Results_annotated[[i]]<-Predict$predict
  #colnames(Results)[i]<-levels(Test.set$peptide)[i]
  #colnames(Results_annotated)[i]<-levels(Test.set$peptide)[i]
  
  Test.set.features[[i]]<-cbind(Test.set.scale, Time=1:(dim(Test.set.scale)[1]))
  Test.set.features[[i]]<-melt(as.data.frame(Test.set.features[[i]]), id.vars = "Time")
  
  g0<-eval(substitute(ggplot(Test.set.features[[i]][-1,], aes(Time, variable)) + 
      geom_tile(aes(fill = value), colour = "white") +
      labs(x = "Time",y = NULL)+
      removeGrid()+
      scale_y_discrete(expand=c(0,0))+
      scale_fill_gradient(low = "white", high = "red",name = "Values")+
      ggtitle(label = levels(Test.set$peptide)[i])+
      theme(legend.position="bottom", panel.background = element_blank(),
            plot.background = element_blank(), plot.margin = unit(c(0.1,0,0,0), "cm"),
            axis.ticks.length = unit(0, "pt"))
      ,list(i = i)))
  interpret.plots[[i]] <- g0
  }
  
  Results<-t(plyr::ldply(Results, rbind))
  colnames(Results)<-levels(Test.set$peptide)
  
  
  Results<-data.frame(RUN=1:(dim(Results)[1]), Results)
  Results_melt <- melt(Results[-1,],id.vars ="RUN")
  decision.map<-ggplot(Results_melt, aes(RUN, variable)) + 
    geom_tile(aes(fill = value), colour = "white") +
    labs(x = "Time",y = NULL)+
    removeGrid()+
    scale_y_discrete(expand=c(0,0))+
    scale_fill_gradient(low = "white", high = "red",name = "Probability")+
    ggtitle(label = "Probability of fail")+
    theme(legend.position="bottom", panel.background = element_blank(),
          plot.background = element_blank(), plot.margin = unit(c(0.1,0,0,0), "cm"),
          axis.ticks.length = unit(0, "pt"))
  
  # gA <- ggplotGrob(g2)
  # gB <- ggplotGrob(g0)
  # maxWidth = grid::unit.pmax(gA$widths[2:5], gB$widths[2:5])
  # gA$widths[2:5] <- as.list(maxWidth)
  # gB$widths[2:5] <- as.list(maxWidth)
  # grid.arrange(gA, ncol=1)  
  
  print(decision.map)
  
  message(paste("Drew the plot for final evaluation"))
  
  print(interpret.plots)
  
  message(paste("Drew the plots for interpretation"))
  
  if (address!=FALSE) {
    dev.off()
  }
}

  
