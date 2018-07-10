simulate_disturbances <- function(guide.set.scale, sim.size){

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
  source("add_features.R")

  for(j in 1:nlevels(guide.set.scale$peptide)){ 
  Data<-c()
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size*5)
  
  Data<-data.frame(idfile=1:(5*sim.size),
                   peptide=rep(levels(guide.set.scale$peptide)[j], (5*sim.size)),
                   sample_data[1], sample_data[2], sample_data[3], sample_data[4])
  Data<- reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("PASS")
  Data <- cbind(Data,RESPONSE)
  
  Data0[[j]]<-Data
  }
Data0<-add_features(guide.set.scale, Data0)

#generate out-of-control observations
#Logarithmic increase 
for(j in 1:nlevels(guide.set$peptide)){
Data<-c()
sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
for(i in 1:sim.size){
  Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]+log(i,base=2)*IQR(sample_data[,1]), 
                       sample_data[i,2]+log(i,base=2)*IQR(sample_data[,2]),
                       sample_data[i,3]+log(i,base=2)*IQR(sample_data[,3]),
                       sample_data[i,4]+log(i,base=2)*IQR(sample_data[,4])))
}
Data<- as.data.frame(Data,stringsAsFactors = F)
for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
RESPONSE<-c(rep("FAIL",sim.size))
Data<- cbind(Data,RESPONSE)
Data1[[j]]<-Data
}
Data1<-add_features(guide.set.scale, Data1)

#Logarithmic decrease 
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]-log(i,base=2)*IQR(sample_data[,1]), 
                       sample_data[i,2]-log(i,base=2)*IQR(sample_data[,2]),
                       sample_data[i,3]-log(i,base=2)*IQR(sample_data[,3]),
                       sample_data[i,4]-log(i,base=2)*IQR(sample_data[,4])))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data2[[j]]<-Data
}
Data2<-add_features(guide.set.scale, Data2)

#generate out-of-control observations for a positive shift 
for(j in 1:nlevels(guide.set.scale$peptide)){
  Data<-c()
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  Data<-data.frame(idfile=1:sim.size,peptide=rep(levels(guide.set.scale$peptide)[j],sim.size),
                   sample_data[1]+2.0*IQR(sample_data[,1]), 
                   sample_data[2]+2.0*IQR(sample_data[,2]), 
                   sample_data[3]+2.0*IQR(sample_data[,3]), 
                   sample_data[4]+2.0*IQR(sample_data[,4]))
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<- reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("FAIL")
  Data <- cbind(Data,RESPONSE)
  Data3[[j]]<-Data
}
Data3<-add_features(guide.set.scale, Data3)

#generate out-of-control observations for a negative shift
for(j in 1:nlevels(guide.set.scale$peptide)){
  Data<-c()
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  Data<-data.frame(idfile=1:sim.size,peptide=rep(levels(guide.set.scale$peptide)[j],sim.size),
                   sample_data[1]-2.0*IQR(sample_data[,1]), 
                   sample_data[2]-2.0*IQR(sample_data[,2]), 
                   sample_data[3]-2.0*IQR(sample_data[,3]), 
                   sample_data[4]-2.0*IQR(sample_data[,4]))
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<- reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("FAIL")
  Data <- cbind(Data,RESPONSE)
  Data4[[j]]<-Data
}
Data4<-add_features(guide.set.scale, Data4)

#Cyclic pattern 

for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]+2*sin(2*pi*i/sim.size), 
                       sample_data[i,2]+2*sin(2*pi*i/sim.size), 
                       sample_data[i,3]+2*sin(2*pi*i/sim.size),
                       sample_data[i,4]+2*sin(2*pi*i/sim.size)))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data5[[j]]<-Data
}
Data5<-add_features(guide.set.scale, Data5)

#Merge all types of disturbances + in-control observations
for(j in 1:nlevels(guide.set.scale$peptide)){
Data.set[[j]]<-rbind(Data0[[j]],Data1[[j]],Data2[[j]],Data3[[j]],Data4[[j]],Data5[[j]])
}
Data.set<-cbind(Data.set[[1]], 
                subset(Data.set[[2]],select = -c(RESPONSE,idfile)), 
                subset(Data.set[[3]],select = -c(RESPONSE,idfile)),
                subset(Data.set[[4]],select = -c(RESPONSE,idfile)),
                subset(Data.set[[5]],select = -c(RESPONSE,idfile)),
                subset(Data.set[[6]],select = -c(RESPONSE,idfile)),
                subset(Data.set[[7]],select = -c(RESPONSE,idfile)),
                subset(Data.set[[8]],select = -c(RESPONSE,idfile)))
return(Data.set)
}