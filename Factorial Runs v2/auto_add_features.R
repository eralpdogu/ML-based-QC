add_features<-function(temp.Data){
    for(i in 1:4){
      v <- numeric(length(temp.Data[,i]))
      MR <- numeric(length(temp.Data[,i]))
      CUSUMpoz.m <- numeric(length(temp.Data[,i]))
      CUSUMneg.m <- numeric(length(temp.Data[,i]))
      CUSUMpoz.v <- numeric(length(temp.Data[,i]))
      CUSUMpoz.m[1]<-0
      CUSUMneg.m[1]<-0
      CUSUMpoz.v[1]<-0
      MR[1]<-0
      v[1] <- (sqrt(abs(temp.Data[1,i]))-0.822)/0.349
      d<-0.5
      for(k in 2:length(temp.Data[,i])) {
        MR[k] <- abs(temp.Data[k,i]-temp.Data[(k-1),i])
        v[k] <- (sqrt(abs(temp.Data[k,i]))-0.822)/0.349
        CUSUMpoz.m[k] <- max(0,(temp.Data[k,i]-d+CUSUMpoz.m[k-1]))
        CUSUMneg.m[k] <- max(0,(-d-temp.Data[k,i]+CUSUMneg.m[k-1]))
        CUSUMpoz.v[k] <- max(0,(v[k]-d+CUSUMpoz.v[k-1]))
      }
      #addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMneg.m, CUSUMpoz.v)
      addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMneg.m)
      colnames(addfeatures) <- paste(colnames(temp.Data)[i], colnames(addfeatures), sep = ".")
      temp.Data<-cbind(temp.Data,addfeatures) # full dataset with external features
    }
  return(temp.Data)
}