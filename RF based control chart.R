setwd("/Users/ed/Dropbox/2. MSstatsQC Paper 2-QCloud")
setwd("C:/Users/aidata-1508/Dropbox/2. MSstatsQC Paper 2-QCloud")

install.packages('caret', dependencies = TRUE) 
install.packages('gtools',dependencies = TRUE)

library(caret)
library(boot)
library(gtools) #smartbin
library(ggplot2)
library(MASS) #mvnrnd
library(pdp)

Train <- read.csv('lumos_training_set.csv')
Test <- read.csv('lumos_all_set.csv')

#remove repeated measurements and reshape the dataset
ind <- which(with( Train, (Train$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | Train$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
S0<-Train[-ind,]
S0<-S0[,-2]
Train<-S0
Train$PepSeq<- gsub("\\(Carbamidomethyl\\)","",Train$PepSeq)
S0$PepSeq<- gsub("\\(Carbamidomethyl\\)","",S0$PepSeq)
S0 <- reshape(S0, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("GO")
S0 <- cbind(S0,RESPONSE)

ind <- which(with( Test, (Test$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | Test$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
Data0<-Test[-ind,]
Data0<-Data0[,-2]
Data0$PepSeq<- gsub("\\(Carbamidomethyl\\)","",Data0$PepSeq)
Data1 <- Data0[1:8 + rep(seq(0, nrow(Data0), by=100), each=8),]
Data1 <- reshape(Data1, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("NOGO")
Data <- cbind(Data1,RESPONSE)
######################################################

N0<-length(S0[,1])
Nw<-10
p0<-c()
p1<-c()
imp<-list()
SW<-c()
Data.model<-c()
Predict<-c()
nTrees<-c()
for (i in Nw:length(Data[,1])) {
  SW<-Data[(i-Nw+1):i,]
  Data.model<-rbind(S0,SW)
  rownames(Data.model)<-1:(N0+Nw)
  Data.model<-Data.model[,-1]
  fitControl <- trainControl(classProbs = T, method = "oob")
  
  fit <- train(as.factor(RESPONSE) ~ ., 
               data = Data.model, 
               method="rf",
               preProcess = c("center", "scale", "nzv"),
               trainControl=fitControl,  
               tuneGrid = data.frame(mtry = 6))
  
  Predict<-predict(fit, type='prob', OOB=TRUE)
  nTrees[i]<-fit$finalModel$ntree
  p1[i]<-sum(Predict[(N0+1):(N0+Nw),2])/Nw
  importance<-varImp(fit)$importance
  imp[[i]]<-t(importance)
  
}

pdpTA <- partial(fit, pred.var = "TotalArea.LVNELTEFAK",plot = TRUE, prob=TRUE)
pdpRT <- partial(fit, pred.var = "RT.LVNELTEFAK", plot = TRUE, prob=TRUE)
pdpMA <- partial(fit, pred.var = "MassAccu.LVNELTEFAK", plot = TRUE, prob=TRUE)
pdpFWHM <- partial(fit, pred.var = "FWHM.LVNELTEFAK", plot = TRUE, prob=TRUE)
grid.arrange(pdpRT, pdpTA, pdpMA, pdpFWHM, ncol = 4)
pdp1<- partial(fit, pred.var = c("TotalArea.LVNELTEFAK","RT.LVNELTEFAK"), plot = TRUE, prob=TRUE, chull=TRUE)
pdp2<- partial(fit, pred.var = c("MassAccu.LVNELTEFAK","RT.LVNELTEFAK"), plot = TRUE, prob=TRUE, chull=TRUE)
pdp3<- partial(fit, pred.var = c("FWHM.LVNELTEFAK", "RT.LVNELTEFAK"), plot = TRUE, prob=TRUE, chull=TRUE)
grid.arrange(pdp1, pdp2, pdp3,ncol = 3)

IMP<-as.data.frame(do.call(rbind, imp))
IMP<-stack(data.frame((sapply(IMP,c))))
IMP<-cbind(IMP, run=rep(Nw:length(Data[,1]),4))
plot2<-ggplot(data=IMP,aes(x=run, y=values, color=ind))+
  geom_line(size=1)+
  geom_point(size=2)+
  ylab("Importance")+
  xlab("run")+
  scale_fill_continuous(guide = guide_legend()) +
  theme(legend.position="bottom")

chart.stat<-as.data.frame(cbind(p1,CL=0.7,run=1:length(p1)))
plot1<-ggplot(data=chart.stat,aes(y=p1, x=run))+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_hline(yintercept=chart.stat$CL, linetype = "dashed", color="red")
cowplot::plot_grid(plot1, plot2, labels = c("RF based control chart", "Importance"), align = "v")

#Simulation 1
#generate multivariate normal data
#parameters from a training sample
n<-100 #incontrol observations 
m<-100 #ooc observations
Data<-c()
Data1<-c()
S0<-c()

#one peptide LVNELTEFAK 
mean <-c(with(data=Train,tapply(RT,INDEX=PepSeq,FUN=mean))[4],
         with(data=Train,tapply(TotalArea,INDEX=PepSeq,FUN=mean)) [4],
         with(data=Train,tapply(MassAccu,INDEX=PepSeq,FUN=mean)) [4],
         with(data=Train,tapply(FWHM,INDEX=PepSeq,FUN=mean)) [4]
         )
covar<-cov(Train[Train$PepSeq=="LVNELTEFAK",c(3,6,7,8)])
S0<-data.frame(idfile=1:n,PepSeq=rep("LVNELTEFAK",n),mvrnorm(n, mean, covar))
colnames(S0)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")
S0 <- reshape(S0, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("GO")
S0 <- cbind(S0,RESPONSE)

Data1<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("LVNELTEFAK",m),u
                  mvrnorm(m, mean+(sqrt(0.0*c(covar[1,1],1.0*covar[2,2],3.0*covar[3,3],1.0*covar[4,4]))), 
                  covar))
colnames(Data1)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")
Data1 <- reshape(Data1, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("NOGO")
Data <- cbind(Data1,RESPONSE)

#Simulation 2
#generate multivariate normal data 4 peptides 
#parameters from a training sample
n<-1000 #incontrol observations 
m<-100 #ooc observations
Data<-c()
Data01<-c()
Data02<-c()
Data03<-c()
Data04<-c()
S01<-c()
S02<-c()
S03<-c()
S04<-c()
S0<-c()

#one peptide EACFAVEGPK 
mean <-c(with(data=Train,tapply(RT,INDEX=PepSeq,FUN=mean))[1],
         with(data=Train,tapply(TotalArea,INDEX=PepSeq,FUN=mean)) [1],
         with(data=Train,tapply(MassAccu,INDEX=PepSeq,FUN=mean)) [1],
         with(data=Train,tapply(FWHM,INDEX=PepSeq,FUN=mean)) [1]
)
covar<-cov(Train[Train$PepSeq=="EACFAVEGPK",c(3,6,7,8)])

S01<-data.frame(idfile=1:n,PepSeq=rep("EACFAVEGPK",n),mvrnorm(n, mean, covar))
colnames(S01)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

Data01<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("EACFAVEGPK",n),
                   mvrnorm(m, mean+(sqrt(0.0*c(covar[1,1],0.0*covar[2,2],0.0*covar[3,3],0.0*covar[4,4]))), 
                           covar))
colnames(Data01)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

#one peptide ECCHGDLLECADDR 
mean <-c(with(data=Train,tapply(RT,INDEX=PepSeq,FUN=mean))[2],
         with(data=Train,tapply(TotalArea,INDEX=PepSeq,FUN=mean)) [2],
         with(data=Train,tapply(MassAccu,INDEX=PepSeq,FUN=mean)) [2],
         with(data=Train,tapply(FWHM,INDEX=PepSeq,FUN=mean)) [2]
)
covar<-cov(Train[Train$PepSeq=="ECCHGDLLECADDR",c(3,6,7,8)])

S02<-data.frame(idfile=1:n,PepSeq=rep("ECCHGDLLECADDR",n), mvrnorm(n, mean, covar))
colnames(S02)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

Data02<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("ECCHGDLLECADDR",n), 
                   mvrnorm(m, mean+(sqrt(0.0*c(covar[1,1],0.0*covar[2,2],0.0*covar[3,3],0.0*covar[4,4]))), 
                           covar))
colnames(Data02)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

#one peptide  HLVDEPQNLIK 
mean <-c(with(data=Train,tapply(RT,INDEX=PepSeq,FUN=mean))[3],
         with(data=Train,tapply(TotalArea,INDEX=PepSeq,FUN=mean)) [3],
         with(data=Train,tapply(MassAccu,INDEX=PepSeq,FUN=mean)) [3],
         with(data=Train,tapply(FWHM,INDEX=PepSeq,FUN=mean)) [3]
)
covar<-cov(Train[Train$PepSeq=="HLVDEPQNLIK",c(3,6,7,8)])

S03<-data.frame(idfile=1:n,PepSeq=rep("HLVDEPQNLIK",n),mvrnorm(n, mean, covar))
colnames(S03)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

Data03<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("HLVDEPQNLIK",n),
                   mvrnorm(m, mean+(sqrt(0.0*c(covar[1,1],0.0*covar[2,2],0.0*covar[3,3],0.0*covar[4,4]))), 
                           covar))
colnames(Data03)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

#one peptide  LVNELTEFAK 
mean <-c(with(data=Train,tapply(RT,INDEX=PepSeq,FUN=mean))[4],
         with(data=Train,tapply(TotalArea,INDEX=PepSeq,FUN=mean)) [4],
         with(data=Train,tapply(MassAccu,INDEX=PepSeq,FUN=mean)) [4],
         with(data=Train,tapply(FWHM,INDEX=PepSeq,FUN=mean)) [4]
)
covar<-cov(Train[Train$PepSeq=="LVNELTEFAK",c(3,6,7,8)])

S04<-data.frame(idfile=1:n,PepSeq=rep("LVNELTEFAK",n),mvrnorm(n, mean, covar))
colnames(S04)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

Data04<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("LVNELTEFAK",n),
                   mvrnorm(m, mean+(sqrt(3.0*c(covar[1,1],3.0*covar[2,2],3.0*covar[3,3],3.0*covar[4,4]))), 
                           covar))
colnames(Data04)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

S0<-rbind(S01,S02,S03,S04)
S0<- reshape(S0, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("GO")
S0 <- cbind(S0,RESPONSE)

Data0<-rbind(Data01,Data02,Data03,Data04)
Data0<-reshape(Data0, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("NOGO")
Data <- cbind(Data0,RESPONSE)