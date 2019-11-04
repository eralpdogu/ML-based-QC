
  Data0<-list()  
  Data1<-list() 
  Data.set<-list()
  #Peptide :- Name of the peptide
  #convert input peptide to character (input without quotes)
  #peptide = as.character(substitute(Peptide))
  #print(peptide)
  
  #generate in-control observations
  source("sample_density_function.R")
  source("auto_add_features.R")
  source("robust_scaling.R")
  
  beta=+4
  sim.size=25
  
  sample_density_sim <- function(guide.set, peptide, n){
    sample_data<-c()
    
    dat.dens = stats::density(guide.set[guide.set$peptide == peptide,3], n=2^10)
    sim.sample.RT = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
    
    dat.dens = stats::density(guide.set[guide.set$peptide == peptide,4], n=2^10)
    sim.sample.TotalArea = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
    
    dat.dens = stats::density(guide.set[guide.set$peptide == peptide,5], n=2^10)
    sim.sample.MassAccu = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
    
    dat.dens = stats::density(guide.set[guide.set$peptide == peptide,6], n=2^10)
    sim.sample.FWHM = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
    
    sample_data <- data.frame(sim.sample.RT,sim.sample.TotalArea,sim.sample.MassAccu, sim.sample.FWHM)
    names(sample_data) <- c("RT", "TotalArea", "MassAccu", "FWHM")
    
    return(sample_data)
  }
  
  for(j in 1:nlevels(guide.set$peptide)){ 
    Data<-c()
    sample_data <- sample_density_sim(guide.set, guide.set$peptide[j], sim.size)
    
    Data<-data.frame(idfile=1:(sim.size),
                     peptide=rep(levels(guide.set$peptide)[j], (sim.size)),
                     sample_data[1], sample_data[2], sample_data[3], sample_data[4])
    RESPONSE<-c("PASS")
    Data <- cbind(Data,RESPONSE)
    Data0[[j]]<-Data
  }
  
  
  #generate out-of-control observations
  #Monotonic increase in RT
  for(j in 1:5){
    Data<-c()
    sample_data <- sample_density_sim(guide.set,guide.set$peptide[j], sim.size)
    
    Data<-data.frame(idfile=(sim.size+1):(sim.size*2),
                     peptide=rep(levels(guide.set$peptide)[j], (sim.size)),
                     sample_data[1], sample_data[2], sample_data[3], sample_data[4])
    RESPONSE<-c("FAIL")
    Data <- cbind(Data,RESPONSE)
    Data1[[j]]<-Data
  }
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  
  for(j in 6:8){
    Data<-c()
    sample_data <- sample_density_sim(guide.set,guide.set$peptide[j], sim.size)
    for(i in 1:sim.size){
      Data<-rbind(Data,c((i+sim.size),rep(levels(guide.set$peptide)[j],1),
                         sample_data[i,1]+ beta*mad(sample_data[,1]),
                         sample_data[i,2],
                         sample_data[i,3],
                         sample_data[i,4]))
    }
    Data<- as.data.frame(Data,stringsAsFactor = F)
    colnames(Data)<-c("idfile", "peptide", colnames(sample_data))
    for (i in c(1,3:ncol(Data))){ Data[,i]<-as.numeric.factor(Data[,i])}
    RESPONSE<-c(rep("FAIL",sim.size))
    Data<- cbind(Data,RESPONSE)
    Data1[[j]]<-Data
  }
  
  #Merge all types of disturbances + in-control observations
  for(j in 1:nlevels(guide.set$peptide)){
    Data.set[[j]]<-rbind(Data0[[j]], Data1[[j]])
  }
  Data.set<-rbind(Data.set[[1]], 
                  Data.set[[2]], 
                  Data.set[[3]],
                  Data.set[[4]],
                  Data.set[[5]],
                  Data.set[[6]],
                  Data.set[[7]],
                  Data.set[[8]])
 
  #Test.set<-Data.set
  
  SimData<-data.frame(Data.set[,1:2], Annotations=NA, Data.set[,3:6])
  colnames(SimData)<-c("Run", "Precursor", "Annotations", "RetentionTime", "TotalArea", "FWHM", "MassAccuracy")
  
  # Sim.set.scale <- SimData[SimData$Precursor==levels(SimData$Precursor)[1],c(4:ncol(SimData))]
  # for(k in 1:ncol(Sim.set.scale)){
  #   Sim.set.scale[,k]=(Sim.set.scale[,k]-median(Sim.set.scale[,k]))/mad(Sim.set.scale[,k])
  # }
  # Sim.set<-Sim.set.scale
  # 
  # for(i in 2:nlevels(SimData$Precursor)){
  #   
  #   Sim.set.scale <- SimData[SimData$Precursor==levels(SimData$Precursor)[i],c(4:ncol(SimData))]
  #   
  #   for(k in 1:ncol(Sim.set.scale)){
  #     Sim.set.scale[,k]=(Sim.set.scale[,k]-median(Sim.set.scale[,k]))/mad(Sim.set.scale[,k])
  #   }
  #   Sim.set<-rbind(Sim.set, Sim.set.scale)
  # }
  # Sim.set<-cbind(SimData[,1:2], Sim.set)
  
  
  #Simdata_melt <- melt(SimData,id.vars =c("Precursor","Run"))
  
  ggplot(SimData, aes(Run, RetentionTime)) + 
    geom_point(size = 0.5)+
    geom_line()+ 
    #geom_smooth(method="lm", col="black")+
    # geom_smooth(data=filter(Simdata_melt, 
    #                         Simdata_melt$variable == "RetentionTime"&Simdata_melt$Run>25), 
    #             aes(Run, value), method = "lm") +
    #geom_point(data=filter(Simdata_melt, Simdata_melt$variable == "RetentionTime"), 
    #           aes(Run, value))+
    #geom_smooth(data=filter(Simdata_melt, Simdata_melt$variable == "RetentionTime"), 
    #            aes(Run, value), method="glm")+
    ylab("Retention Time")+
    xlab("Time")+
    facet_wrap(~Precursor,scales = "free", ncol = 4)+
    scale_color_manual(values = c("#F0E442", "#0072B2", "#CC79A7", "#D55E00"))+
    labs(color = "Metric")+
    theme(legend.position="bottom", panel.background = element_blank(),
          plot.background = element_blank(), plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),
          axis.ticks.length = unit(0, "pt"))
  
  
  MSstatsQC.ML.testR(Data.set[,1:6], guide.set)
  
  write.csv(Data.set, "SimData.csv")
