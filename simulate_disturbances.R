simulate_disturbances <- function(guide.set.scale, sim.size){

n<-sim.size#incontrol observations 
Data.set<-c()
Data0<-list() #storage for in-control observations 
Data1<-list() #storage for logarithmic increase 
Data2<-list() #storage for logarithmic decrease
Data3<-list() #storage for mean increase
Data4<-list() #storage for mean decreaase
Data5<-list() #storage for cyclic pattern
  
Data.set<-list()

  #Peptide :- Name of the peptide
  #convert input peptide to character (input without quotes)
  #peptide = as.character(substitute(Peptide))
  #print(peptide)
  
  #generate in-control observations
  source("sample_density_function.R")

  for(j in 1:nlevels(guide.set.scale$peptide)){

  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], n)
  
  Data<-data.frame(idfile=1:n,peptide=rep(levels(guide.set.scale$peptide)[j], n),
                 sample_data$RT, sample_data$TotalArea, sample_data$MassAccu, sample_data$FWHM)
  colnames(Data)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  Data<- reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("PASS")
  Data <- cbind(Data,RESPONSE)
  
  Data0[[j]]<-Data
  }

for(j in 1:nlevels(guide.set.scale$peptide)){
  for(i in 2:5){
    v <- numeric(length(Data0[[j]][,i]))
    MR <- numeric(length(Data0[[j]][,i]))
    CUSUMpoz.m <- numeric(length(Data0[[j]][,i]))
    CUSUMneg.m <- numeric(length(Data0[[j]][,i]))
    CUSUMpoz.v <- numeric(length(Data0[[j]][,i]))
    CUSUMpoz.m[1]<-0
    CUSUMneg.m[1]<-0
    CUSUMpoz.v[1]<-0
    MR[1]<-0
    v[1]<-0
    d<-0.5
    for(k in 2:length(Data0[[j]][,i])) {
      MR[k] <- abs(Data0[[j]][k,i]-Data0[[j]][(k-1),i])
      v[k] <- (sqrt(abs(Data0[[j]][k,i]))-0.822)/0.349
      CUSUMpoz.m[k] <- max(0,(Data0[[j]][k,i]-(d)+CUSUMpoz.m[k-1]))
      CUSUMneg.m[k] <- max(0,((-d)-Data0[[j]][k,i]+CUSUMneg.m[k-1]))
      CUSUMpoz.v[k] <- max(0,(v[k]-(k)+CUSUMpoz.v[k-1]))
    }
    addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMpoz.v)
    colnames(addfeatures) <- paste(levels(guide.set$peptide)[j], colnames(addfeatures), sep = ".")
    Data0[[j]]<-cbind(Data0[[j]],addfeatures) # full dataset with external features
  }
}

#generate out-of-control observations
#Logarithmic increase in FWHM 
for(j in 1:nlevels(guide.set$peptide)){
Data<-c()
sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], n)
for(i in 1:n){
  Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data$RT[i], sample_data$TotalArea[i], sample_data$MassAccu[i],
                       sample_data$FWHM[i]+log(i,base=2)*IQR(sample_data$FWHM)))
}
colnames(Data)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
Data<- as.data.frame(Data,stringsAsFactors = F)
for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
RESPONSE<-c(rep("FAIL",n))
Data<- cbind(Data,RESPONSE)
Data1[[j]]<-Data
}

for(j in 1:nlevels(guide.set$peptide)){
  for(i in 2:5){
  v <- numeric(length(Data1[[j]][,i]))
  MR <- numeric(length(Data1[[j]][,i]))
  CUSUMpoz.m <- numeric(length(Data1[[j]][,i]))
  CUSUMneg.m <- numeric(length(Data1[[j]][,i]))
  CUSUMpoz.v <- numeric(length(Data1[[j]][,i]))
  CUSUMpoz.m[1]<-0
  CUSUMneg.m[1]<-0
  CUSUMpoz.v[1]<-0
  MR[1]<-0
  v[1]<-0
  d<-0.5
  for(k in 2:length(Data1[[j]][,i])) {
    MR[k] <- abs(Data1[[j]][k,i]-Data1[[j]][(k-1),i])
    v[k] <- (sqrt(abs(Data1[[j]][k,i]))-0.822)/0.349
    CUSUMpoz.m[k] <- max(0,(Data1[[j]][k,i]-(d)+CUSUMpoz.m[k-1]))
    CUSUMneg.m[k] <- max(0,((-d)-Data1[[j]][k,i]+CUSUMneg.m[k-1]))
    CUSUMpoz.v[k] <- max(0,(v[k]-(k)+CUSUMpoz.v[k-1]))
  }
  addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMpoz.v)
  colnames(addfeatures) <- paste(levels(guide.set$peptide)[j], colnames(addfeatures), sep = ".")
  Data1[[j]]<-cbind(Data1[[j]],addfeatures) # full dataset with external features
  }
  }

#Logarithmic decrease in FWHM 
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], n)
  for(i in 1:n){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data$RT[i], sample_data$TotalArea[i], sample_data$MassAccu[i],
                       sample_data$FWHM[i]-log(i,base=2)*IQR(sample_data$FWHM)))
  }
  colnames(Data)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",n))
  Data<- cbind(Data,RESPONSE)
  Data2[[j]]<-Data
}

for(j in 1:nlevels(guide.set.scale$peptide)){
  for(i in 2:5){
    v <- numeric(length(Data2[[j]][,i]))
    MR <- numeric(length(Data2[[j]][,i]))
    CUSUMpoz.m <- numeric(length(Data2[[j]][,i]))
    CUSUMneg.m <- numeric(length(Data2[[j]][,i]))
    CUSUMpoz.v <- numeric(length(Data2[[j]][,i]))
    CUSUMpoz.m[1]<-0
    CUSUMneg.m[1]<-0
    CUSUMpoz.v[1]<-0
    MR[1]<-0
    v[1]<-0
    d<-0.5
    for(k in 2:length(Data2[[j]][,i])) {
      MR[k] <- abs(Data2[[j]][k,i]-Data0[[j]][(k-1),i])
      v[k] <- (sqrt(abs(Data2[[j]][k,i]))-0.822)/0.349
      CUSUMpoz.m[k] <- max(0,(Data2[[j]][k,i]-(d)+CUSUMpoz.m[k-1]))
      CUSUMneg.m[k] <- max(0,((-d)-Data2[[j]][k,i]+CUSUMneg.m[k-1]))
      CUSUMpoz.v[k] <- max(0,(v[k]-(k)+CUSUMpoz.v[k-1]))
    }
    addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMpoz.v)
    colnames(addfeatures) <- paste(levels(guide.set$peptide)[j], colnames(addfeatures), sep = ".")
    Data2[[j]]<-cbind(Data2[[j]],addfeatures) # full dataset with external features
  }
}

#generate out-of-control observations for a +2 IQR shift in Total.area---large shift
for(j in 1:nlevels(guide.set.scale$peptide)){
  
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], n)
  
  Data<-data.frame(idfile=1:n,peptide=rep(levels(guide.set.scale$peptide)[j],n),
                   sample_data$RT, sample_data$TotalArea+2.0*IQR(sample_data$TotalArea), 
                   sample_data$MassAccu, sample_data$FWHM)
  colnames(Data)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  Data<- reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("FAIL")
  Data <- cbind(Data,RESPONSE)
  Data3[[j]]<-Data
}

for(j in 1:nlevels(guide.set.scale$peptide)){
  for(i in 2:5){
    v <- numeric(length(Data3[[j]][,i]))
    MR <- numeric(length(Data3[[j]][,i]))
    CUSUMpoz.m <- numeric(length(Data3[[j]][,i]))
    CUSUMneg.m <- numeric(length(Data3[[j]][,i]))
    CUSUMpoz.v <- numeric(length(Data3[[j]][,i]))
    CUSUMpoz.m[1]<-0
    CUSUMneg.m[1]<-0
    CUSUMpoz.v[1]<-0
    MR[1]<-0
    v[1]<-0
    d<-0.5
    for(k in 2:length(Data3[[j]][,i])) {
      MR[k] <- abs(Data3[[j]][k,i]-Data3[[j]][(k-1),i])
      v[k] <- (sqrt(abs(Data3[[j]][k,i]))-0.822)/0.349
      CUSUMpoz.m[k] <- max(0,(Data3[[j]][k,i]-(d)+CUSUMpoz.m[k-1]))
      CUSUMneg.m[k] <- max(0,((-d)-Data3[[j]][k,i]+CUSUMneg.m[k-1]))
      CUSUMpoz.v[k] <- max(0,(v[k]-(k)+CUSUMpoz.v[k-1]))
    }
    addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMpoz.v)
    colnames(addfeatures) <- paste(levels(guide.set$peptide)[j], colnames(addfeatures), sep = ".")
    Data3[[j]]<-cbind(Data3[[j]],addfeatures) # full dataset with external features
  }
}

#generate out-of-control observations for a -2 IQR shift in Total.area---large shift
for(j in 1:nlevels(guide.set.scale$peptide)){
  
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], n)
  
  Data<-data.frame(idfile=1:n,peptide=rep(levels(guide.set.scale$peptide)[j],n),
                   sample_data$RT, sample_data$TotalArea+2.0*IQR(sample_data$TotalArea), 
                   sample_data$MassAccu, sample_data$FWHM)
  colnames(Data)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  Data<- reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("FAIL")
  Data <- cbind(Data,RESPONSE)
  Data4[[j]]<-Data
}

for(j in 1:nlevels(guide.set.scale$peptide)){
  for(i in 2:5){
    v <- numeric(length(Data4[[j]][,i]))
    MR <- numeric(length(Data4[[j]][,i]))
    CUSUMpoz.m <- numeric(length(Data4[[j]][,i]))
    CUSUMneg.m <- numeric(length(Data4[[j]][,i]))
    CUSUMpoz.v <- numeric(length(Data4[[j]][,i]))
    CUSUMpoz.m[1]<-0
    CUSUMneg.m[1]<-0
    CUSUMpoz.v[1]<-0
    MR[1]<-0
    v[1]<-0
    d<-0.5
    for(k in 2:length(Data4[[j]][,i])) {
      MR[k] <- abs(Data4[[j]][k,i]-Data4[[j]][(k-1),i])
      v[k] <- (sqrt(abs(Data4[[j]][k,i]))-0.822)/0.349
      CUSUMpoz.m[k] <- max(0,(Data4[[j]][k,i]-(d)+CUSUMpoz.m[k-1]))
      CUSUMneg.m[k] <- max(0,((-d)-Data4[[j]][k,i]+CUSUMneg.m[k-1]))
      CUSUMpoz.v[k] <- max(0,(v[k]-(k)+CUSUMpoz.v[k-1]))
    }
    addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMpoz.v)
    colnames(addfeatures) <- paste(levels(guide.set$peptide)[j], colnames(addfeatures), sep = ".")
    Data4[[j]]<-cbind(Data4[[j]],addfeatures) # full dataset with external features
  }
}

#Cyclic pattern in MassAccu

for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], n)
  for(i in 1:n){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data$RT[i], sample_data$TotalArea[i], 
                       sample_data$MassAccu[i]+2*sin(2*pi*i/n),
                       sample_data$FWHM[i]))
  }
  colnames(Data)<-c("idfile","peptide","RT","TotalArea","MassAccu","FWHM")
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",n))
  Data<- cbind(Data,RESPONSE)
  Data5[[j]]<-Data
}

for(j in 1:nlevels(guide.set.scale$peptide)){
  for(i in 2:5){
    v <- numeric(length(Data5[[j]][,i]))
    MR <- numeric(length(Data5[[j]][,i]))
    CUSUMpoz.m <- numeric(length(Data5[[j]][,i]))
    CUSUMneg.m <- numeric(length(Data5[[j]][,i]))
    CUSUMpoz.v <- numeric(length(Data5[[j]][,i]))
    CUSUMpoz.m[1]<-0
    CUSUMneg.m[1]<-0
    CUSUMpoz.v[1]<-0
    MR[1]<-0
    v[1]<-0
    d<-0.5
    for(k in 2:length(Data5[[j]][,i])) {
      MR[k] <- abs(Data5[[j]][k,i]-Data5[[j]][(k-1),i])
      v[k] <- (sqrt(abs(Data5[[j]][k,i]))-0.822)/0.349
      CUSUMpoz.m[k] <- max(0,(Data5[[j]][k,i]-(d)+CUSUMpoz.m[k-1]))
      CUSUMneg.m[k] <- max(0,((-d)-Data5[[j]][k,i]+CUSUMneg.m[k-1]))
      CUSUMpoz.v[k] <- max(0,(v[k]-(k)+CUSUMpoz.v[k-1]))
    }
    addfeatures<-cbind(MR,CUSUMpoz.m,CUSUMpoz.v)
    colnames(addfeatures) <- paste(levels(guide.set$peptide)[j], colnames(addfeatures), sep = ".")
    Data5[[j]]<-cbind(Data5[[j]],addfeatures) # full dataset with external features
  }
}

#Merge all types of disturbances + in-control observations
for(j in 1:nlevels(guide.set.scale$peptide)){
Data.set[[j]]<-rbind(Data0[[j]],Data1[[j]],Data2[[j]],Data3[[j]],Data4[[j]],Data5[[j]])
}
return(Data.set)
}