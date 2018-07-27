add_features<-function(guide.set.scale, temp.Data, pep.index){
  for(j in pep.index){
  for(i in 2:5){
    v <- numeric(length(temp.Data[[j]][,i]))
    MR <- numeric(length(temp.Data[[j]][,i]))
    CUSUMpoz.m <- numeric(length(temp.Data[[j]][,i]))
    CUSUMneg.m <- numeric(length(temp.Data[[j]][,i]))
    CUSUMpoz.v <- numeric(length(temp.Data[[j]][,i]))
    CUSUMpoz.m[1]<-0
    CUSUMneg.m[1]<-0
    CUSUMpoz.v[1]<-0
    MR[1]<-0
    v[1]<-0
    d<-0.5
    for(k in 2:length(temp.Data[[j]][,i])) {
      MR[k] <- abs(temp.Data[[j]][k,i]-temp.Data[[j]][(k-1),i])
      v[k] <- (sqrt(abs(temp.Data[[j]][k,i]))-0.822)/0.349
      CUSUMpoz.m[k] <- max(0,(temp.Data[[j]][k,i]-(d)+CUSUMpoz.m[k-1]))
      CUSUMneg.m[k] <- max(0,(-d-temp.Data[[j]][k,i]+CUSUMneg.m[k-1]))
      CUSUMpoz.v[k] <- max(0,(v[k]-d+CUSUMpoz.v[k-1]))
    }
    addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMneg.m, CUSUMpoz.v)
    colnames(addfeatures) <- paste(levels(guide.set.scale$peptide)[j],colnames(temp.Data[[j]])[i], colnames(addfeatures), sep = ".")
    temp.Data[[j]]<-cbind(temp.Data[[j]],addfeatures) # full dataset with external features
  }
  }
  return(temp.Data)
}