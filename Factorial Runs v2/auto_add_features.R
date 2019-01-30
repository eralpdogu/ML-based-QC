add_features<-function(temp.Data){
    for(i in 1:ncol(temp.Data)){
      v <- numeric(length(temp.Data[,i]))
      Moving.range <- numeric(length(temp.Data[,i]))
      Mean.increase <- numeric(length(temp.Data[,i]))
      Mean.decrease <- numeric(length(temp.Data[,i]))
      Var.increase <- numeric(length(temp.Data[,i]))
      Mean.increase[1]<-0
      Mean.decrease[1]<-0
      Var.increase[1]<-0
      Moving.range[1]<-0
      v[1] <- (sqrt(abs(temp.Data[1,i]))-0.822)/0.349
      d<-0.5
      for(k in 2:length(temp.Data[,i])) {
        Moving.range[k] <- abs(temp.Data[k,i]-temp.Data[(k-1),i])
        v[k] <- (sqrt(abs(temp.Data[k,i]))-0.822)/0.349
        Mean.increase[k] <- max(0,(temp.Data[k,i]-d+Mean.increase[k-1]))
        Mean.decrease[k] <- max(0,(-d-temp.Data[k,i]+Mean.decrease[k-1]))
        Var.increase[k] <- max(0,(v[k]-d+Var.increase[k-1]))
      }
      addfeatures<-cbind(Moving.range,Mean.increase,Mean.decrease, Var.increase)
      #addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMneg.m)
      colnames(addfeatures) <- paste(colnames(temp.Data)[i], colnames(addfeatures), sep = ".")
      temp.Data<-cbind(temp.Data,addfeatures) # full dataset with external features
    }
  return(temp.Data)
}
