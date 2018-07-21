simulate_disturbances <- function(guide.set.scale, sim.size){

Data0<-list()  
Data1<-list() 
Data2<-list() 
Data3<-list() 
Data4<-list() 
Data5<-list() 
Data6<-list() 
Data7<-list() 
Data8<-list() 
Data9<-list() 
Data10<-list()
Data11<-list()
Data12<-list() 
Data13<-list() 
Data14<-list() 
Data15<-list()
Data16<-list()
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
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size*16)
  
  Data<-data.frame(idfile=1:(16*sim.size),
                   peptide=rep(levels(guide.set.scale$peptide)[j], (16*sim.size)),
                   sample_data[1], sample_data[2], sample_data[3], sample_data[4])
  Data<- reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c("PASS")
  Data <- cbind(Data,RESPONSE)
  
  Data0[[j]]<-Data
  }
Data0<-add_features(guide.set.scale, Data0)

#generate out-of-control observations
#Monotonic increase in RT
for(j in 1:nlevels(guide.set$peptide)){
Data<-c()
beta=runif(1,0,3)
sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
for(i in 1:sim.size){
  Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]+beta*(i-sim.size), 
                       sample_data[i,2],
                       sample_data[i,3],
                       sample_data[i,4]))
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

#Monotonic increase in Total Area
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2]+beta*(i-sim.size),
                       sample_data[i,3],
                       sample_data[i,4]))
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

#Monotonic increase in Mass Accu
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2],
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data3[[j]]<-Data
}
Data3<-add_features(guide.set.scale, Data3)

#generate out-of-control observations for a negative shift
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2],
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data4[[j]]<-Data
}
Data4<-add_features(guide.set.scale, Data4)

#NO CHANGE Run:1

for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2], 
                       sample_data[i,3],
                       sample_data[i,4]))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("PASS",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data5[[j]]<-Data
}
Data5<-add_features(guide.set.scale, Data5)

#RUN 5: AB
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]+beta*(i-sim.size), 
                       sample_data[i,2]+beta*(i-sim.size), 
                       sample_data[i,3],
                       sample_data[i,4]))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data6[[j]]<-Data
}
Data6<-add_features(guide.set.scale, Data6)

#RUN 6: AC
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]+beta*(i-sim.size), 
                       sample_data[i,2], 
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data7[[j]]<-Data
}
Data7<-add_features(guide.set.scale, Data7)

#RUN 7: BC
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2]+beta*(i-sim.size), 
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data8[[j]]<-Data
}
Data8<-add_features(guide.set.scale, Data8)

#RUN 8: ABC
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]+beta*(i-sim.size), 
                       sample_data[i,2]+beta*(i-sim.size), 
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data9[[j]]<-Data
}
Data9<-add_features(guide.set.scale, Data9)

#RUN 10: AD
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]-beta*(i-sim.size), 
                       sample_data[i,2], 
                       sample_data[i,3],
                       sample_data[i,4]-beta*(i-sim.size)))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data10[[j]]<-Data
}
Data10<-add_features(guide.set.scale, Data10)

#RUN 11: BD
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2]+beta*(i-sim.size), 
                       sample_data[i,3],
                       sample_data[i,4]+beta*(i-sim.size)))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data11[[j]]<-Data
}
Data11<-add_features(guide.set.scale, Data11)

#RUN 12: ABD
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]+beta*(i-sim.size), 
                       sample_data[i,2]+beta*(i-sim.size), 
                       sample_data[i,3],
                       sample_data[i,4]+beta*(i-sim.size)))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data12[[j]]<-Data
}
Data12<-add_features(guide.set.scale, Data12)

#RUN 13: CD
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2], 
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]+beta*(i-sim.size)))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data13[[j]]<-Data
}
Data13<-add_features(guide.set.scale, Data13)

#RUN 14: ACD
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1]+beta*(i-sim.size), 
                       sample_data[i,2], 
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]-beta*(i-sim.size)))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data14[[j]]<-Data
}
Data14<-add_features(guide.set.scale, Data14)

#RUN 15: BCD
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2]+beta*(i-sim.size), 
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]+beta*(i-sim.size)))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data15[[j]]<-Data
}
Data15<-add_features(guide.set.scale, Data15)

#RUN 16: ABCD
for(j in 1:nlevels(guide.set$peptide)){
  Data<-c()
  beta=runif(1,0,3)
  sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
  for(i in 1:sim.size){
    Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                       sample_data[i,1], 
                       sample_data[i,2]+beta*(i-sim.size), 
                       sample_data[i,3]+beta*(i-sim.size),
                       sample_data[i,4]+beta*(i-sim.size)))
  }
  Data<- as.data.frame(Data,stringsAsFactors = F)
  for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric(Data[,i])}
  colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
  Data<-reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
  RESPONSE<-c(rep("FAIL",sim.size))
  Data<- cbind(Data,RESPONSE)
  Data16[[j]]<-Data
}
Data16<-add_features(guide.set.scale, Data16)

#Merge all types of disturbances + in-control observations
for(j in 1:nlevels(guide.set.scale$peptide)){
Data.set[[j]]<-rbind(Data0[[j]],Data1[[j]],Data2[[j]],Data3[[j]],Data4[[j]],Data5[[j]],
                     Data6[[j]],Data7[[j]],Data8[[j]],Data9[[j]],Data10[[j]],Data11[[j]],
                     Data13[[j]],Data14[[j]],Data15[[j]],Data16[[j]])
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
