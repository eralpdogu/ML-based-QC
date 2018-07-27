simulate_drifts_all <-function(guide.set.scale, sim.size){
  
  ####DISTURBANCES for ALL PEPTIDES################################################
  #################################################################################
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
  
  num.pep<-nlevels(guide.set.scale$peptide)
  
  for(j in 1:num.pep){
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

  Data0<-add_features(guide.set.scale, Data0, 1:num.pep)

  for(j in 1:num.pep){Data.set[[j]]<-rbind(Data0[[j]])
  levels(Data.set[[j]][,6])<-c("PASS","FAIL")}

  #generate out-of-control observations
  #Monotonic increase in RT
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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

  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(16*sim.size+1):(16*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(16*sim.size+1):(16*sim.size+sim.size),]<-Data1[[j]]}

  #Monotonic increase in Total Area
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(17*sim.size+1):(17*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(17*sim.size+1):(17*sim.size+sim.size),]<-Data1[[j]]}

  #Monotonic increase in Mass Accu
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }

  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(18*sim.size+1):(18*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(18*sim.size+1):(18*sim.size+sim.size),]<-Data1[[j]]}

  #generate out-of-control observations for a negative shift
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(19*sim.size+1):(19*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(19*sim.size+1):(19*sim.size+sim.size),]<-Data1[[j]]}

  #NO CHANGE Run:1

  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(20*sim.size+1):(20*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(20*sim.size+1):(20*sim.size+sim.size),]<-Data1[[j]]}

  #RUN 5: AB
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(21*sim.size+1):(21*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(21*sim.size+1):(21*sim.size+sim.size),]<-Data1[[j]]}

  #RUN 6: AC
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(22*sim.size+1):(22*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(22*sim.size+1):(22*sim.size+sim.size),]<-Data1[[j]]}


  #RUN 7: BC
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(23*sim.size+1):(23*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(23*sim.size+1):(23*sim.size+sim.size),]<-Data1[[j]]}

  #RUN 8: ABC
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(24*sim.size+1):(24*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(24*sim.size+1):(24*sim.size+sim.size),]<-Data1[[j]]}

  #RUN 10: AD
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(25*sim.size+1):(25*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(25*sim.size+1):(25*sim.size+sim.size),]<-Data1[[j]]}

  #RUN 11: BD
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(26*sim.size+1):(26*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(26*sim.size+1):(26*sim.size+sim.size),]<-Data1[[j]]}


  #RUN 12: ABD
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(27*sim.size+1):(27*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(27*sim.size+1):(27*sim.size+sim.size),]<-Data1[[j]]}


  #RUN 13: CD
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(28*sim.size+1):(28*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(28*sim.size+1):(28*sim.size+sim.size),]<-Data1[[j]]}


  #RUN 14: ACD
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(29*sim.size+1):(29*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(29*sim.size+1):(29*sim.size+sim.size),]<-Data1[[j]]}

  #RUN 15: BCD
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(30*sim.size+1):(30*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(30*sim.size+1):(30*sim.size+sim.size),]<-Data1[[j]]}


  #RUN 16: ABCD
  for(j in 1:num.pep){
    Data<-c()
    sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[j], sim.size)
    for(i in 1:sim.size){
      beta=runif(1,0,3)
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
    Data1[[j]]<-Data
  }
  Data1<-add_features(guide.set.scale, Data1, 1:num.pep)
  Data1[[j]][,1]<-(31*sim.size+1):(31*sim.size+sim.size)

  for(j in 1:num.pep){Data.set[[j]][(31*sim.size+1):(31*sim.size+sim.size),]<-Data1[[j]]}
  Data.set.temp<-cbind(Data.set[[1]],
                  subset(Data.set[[2]],select = -c(RESPONSE,idfile)))

  for(j in 3:num.pep){Data.set.temp<-cbind(Data.set.temp,
                                subset(Data.set[[j]],select = -c(RESPONSE,idfile)))}
  Data.set<-Data.set.temp

}