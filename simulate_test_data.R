simulate_test_data <- function(guide.set.scale, sim.size, beta){
  
  Data0<-list()  
  Data1<-list() 
  Data.set<-list()
  #Peptide :- Name of the peptide
  #convert input peptide to character (input without quotes)
  #peptide = as.character(substitute(Peptide))
  #print(peptide)
  
  #generate in-control observations
  source("sample_density_function.R")
  source("add_features.R")
  set.seed(123)
  for(j in 1:nlevels(guide.set.scale$peptide)){ 
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    
    Data<-data.frame(idfile=1:(sim.size),
                     peptide=rep(levels(guide.set.scale$peptide)[j], (sim.size)),
                     sample_data[1], sample_data[2], sample_data[3], sample_data[4])
    Data<- reshape(Data, idvar = "idfile", timevar = "peptide", direction = "wide")
    RESPONSE<-c("PASS")
    Data <- cbind(Data,RESPONSE)
    Data0[[j]]<-Data
  }
  Data0<-add_features(guide.set.scale, Data0, pep.index = 1:num.pep)
  
  #generate out-of-control observations
  #Monotonic increase in RT
  for(j in 1:nlevels(guide.set$peptide)){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      Data<-rbind(Data,c(i,rep(levels(guide.set.scale$peptide)[j],1),
                         sample_data[i,1]-beta*(i-sim.size),
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
  Data1<-add_features(guide.set.scale, Data1, pep.index = 1:num.pep)
  
  #Merge all types of disturbances + in-control observations
  for(j in 1:nlevels(guide.set.scale$peptide)){
    Data.set[[j]]<-rbind(Data0[[j]],Data1[[j]])
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