setwd("/Users/ed/Dropbox/2. MSstatsQC Paper 2-QCloud")
setwd("C:/Users/aidata-1508/Dropbox/2. MSstatsQC Paper 2-QCloud")

install.packages('caret', dependencies = TRUE) 
install.packages('gtools',dependencies = TRUE)
library(caret)
library(boot)
library(gtools) #smartbin
library(reshape2) #melt
library(ggplot2)

Train <- read.csv('lumos_training_set.csv')
Test <- read.csv('lumos_all_set.csv')

#remove repeated measurements and reshape the dataset
ind <- which(with( Train, (Train$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | Train$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
S0<-Train[-ind,]
S0<-S0[,-2]
Train<-S0
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
Nw<-2
p0<-c()
p1<-c()
imp<-list()

for (i in N0+1:length(S0[,1])){
  SW<-Data[(i-N0):(i+Nw-1-N0),]
  Data.model<-rbind(S0,SW)
  Data.model<-Data.model[,-1]
  
  fit <- train(as.factor(RESPONSE) ~ ., data = Data.model, method="rf",  preProcess = c("center", "scale", "nzv")
               , trainControl=trainControl( method="oob" ))
  
  Predict<-predict(fit, type='prob')
  p0[i]<-sum(Predict[,1])/length(S0[,1])
  importance<-varImp(fit)$importance
  imp[[i]]<-t(importance)

}

IMP<-as.data.frame(melt(do.call(smartbind,imp)))
IMP<-data.frame(c(IMP,rep(1:i,length(importance$Overall))))

plot(p0, type = 'l')
abline(h = CL,col = "gray60")

ggplot(data=IMP,aes(y=value, x=1:340, col=variable))+geom_line

#Simulation 1
#generate multivariate normal data
#parameters from a training sample
n<-10 #incontrol observations 
m<-10 #ooc observations
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

Data1<-data.frame(idfile=(1+m):(n+m),PepSeq=rep("LVNELTEFAK",m),mvrnorm(m, mean+3*c(covar[1],covar[2],covar[3],covar[4]), covar))
colnames(Data1)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")
Data1 <- reshape(Data1, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("NOGO")
Data <- cbind(Data1,RESPONSE)

#Simulation 2
#generate multivariate normal data 4 peptides 
#parameters from a training sample
n<-10 #incontrol observations 
m<-10 #ooc observations
#one peptide EACFAVEGPK 
mean <-c(with(data=Train,tapply(RT,INDEX=PepSeq,FUN=mean))[1],
         with(data=Train,tapply(TotalArea,INDEX=PepSeq,FUN=mean)) [1],
         with(data=Train,tapply(MassAccu,INDEX=PepSeq,FUN=mean)) [1],
         with(data=Train,tapply(FWHM,INDEX=PepSeq,FUN=mean)) [1]
)
covar<-cov(Train[Train$PepSeq=="EACFAVEGPK",c(3,6,7,8)])

S01<-data.frame(idfile=1:n,PepSeq=rep("EACFAVEGPK",n),mvrnorm(n, mean, covar))
colnames(S01)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

Data01<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("EACFAVEGPK",n),mvrnorm(m, mean+3*c(covar[1],covar[2],covar[3],covar[4]), covar))
colnames(Data01)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

#one peptide ECCHGDLLECADDR 
mean <-c(with(data=Train,tapply(RT,INDEX=PepSeq,FUN=mean))[2],
         with(data=Train,tapply(TotalArea,INDEX=PepSeq,FUN=mean)) [2],
         with(data=Train,tapply(MassAccu,INDEX=PepSeq,FUN=mean)) [2],
         with(data=Train,tapply(FWHM,INDEX=PepSeq,FUN=mean)) [2]
)
covar<-cov(Train[Train$PepSeq=="ECCHGDLLECADDR",c(3,6,7,8)])

S02<-data.frame(idfile=1:n,PepSeq=rep("ECCHGDLLECADDR",n),mvrnorm(n, mean, covar))
colnames(S02)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

Data02<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("ECCHGDLLECADDR",n),mvrnorm(m, mean+3*c(covar[1],covar[2],covar[3],covar[4]), covar))
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

Data03<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("HLVDEPQNLIK",n),mvrnorm(m, mean+3*c(covar[1],covar[2],covar[3],covar[4]), covar))
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

Data04<-data.frame(idfile=(n+1):(n+m),PepSeq=rep("LVNELTEFAK",n),mvrnorm(m, mean+3*c(covar[1],covar[2],covar[3],covar[4]), covar))
colnames(Data04)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

S0<-rbind(S01,S02,S03,S04)
S0<- reshape(S0, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("GO")
S0 <- cbind(S0,RESPONSE)

Data0<-rbind(Data01,Data02,Data03,Data04)
Data0<-reshape(Data0, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("NOGO")
Data <- cbind(Data0,RESPONSE)