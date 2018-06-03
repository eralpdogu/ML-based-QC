setwd("/Users/ed/Dropbox/2. MSstatsQC Paper 2-QCloud")
setwd("C:/Users/aidata-1508/Dropbox/2. MSstatsQC Paper 2-QCloud")

install.packages('caret', dependencies = TRUE) 
install.packages('gtools',dependencies = TRUE)
install.packages('pdp',dependencies = TRUE)

library(caret)
library(boot)
library(gtools) #smartbin
library(ggplot2)
library(MASS) #mvnrnd
library(pdp)
library(dplyr) #sample_n

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

#Simulation 1
#generate multivariate normal data
#parameters from a training sample
n<-10000 #incontrol observations 
Data<-c()
Data0<-c()
Data1<-c()
Data2<-c()
Data3<-c()
Data4<-c()
S0<-c()

#one peptide LVNELTEFAK 
mean <-c(with(data=Train,tapply(RT,INDEX=PepSeq,FUN=mean))[4],
         with(data=Train,tapply(TotalArea,INDEX=PepSeq,FUN=mean)) [4],
         with(data=Train,tapply(MassAccu,INDEX=PepSeq,FUN=mean)) [4],
         with(data=Train,tapply(FWHM,INDEX=PepSeq,FUN=mean)) [4]
)
covar<-cov(Train[Train$PepSeq=="LVNELTEFAK",c(3,6,7,8)])
#generate in-control observations

S0<-data.frame(idfile=1:(4*n),PepSeq=rep("LVNELTEFAK",n),mvrnorm(n, mean, covar))
colnames(S0)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")
S0<- reshape(S0, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("GO")
S0 <- cbind(S0,RESPONSE)
#Logarithmic drift
Data1<-data.frame(idfile=((1):(n)),PepSeq=rep("LVNELTEFAK",n),
                  mvrnorm(n, mean+c(1.0*sqrt(covar[1,1]),3.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),1.0*sqrt(covar[4,4])), 
                          covar))
# Data1<-as.data.frame(matrix(0,n,6))
# for(i in 1:n){
#   Data1[i,]<-c(idfile=i,PepSeq=rep("LVNELTEFAK",1),
#                mvrnorm(1, mean+c(1.0*sqrt(covar[1,1]),1.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),1.0*sqrt(covar[4,4])), 
#                        covar))
# }
colnames(Data1)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

#generate out-of-control observations for a 3 sigma step shift in Total.area---large shift
Data2<-data.frame(idfile=((n+1):(2*n)),PepSeq=rep("LVNELTEFAK",n),
                  mvrnorm(n, mean+c(1.0*sqrt(covar[1,1]),3.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),1.0*sqrt(covar[4,4])), 
                          covar))
colnames(Data2)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

#generate out-of-control observations for a 2 sigma fluctuation in Mass.accu---large shift
covar1<-covar
covar1[3,3]<-3*covar[3,3]
Data3<-data.frame(idfile=((2*n+1):(3*n)),PepSeq=rep("LVNELTEFAK",n),
                  mvrnorm(n, mean, 
                          covar1))
colnames(Data3)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

#generate out-of-control observations for a 1.5 sigma step shift in RT---small shift
Data4<-data.frame(idfile=(3*n+1):(4*n),PepSeq=rep("LVNELTEFAK",n),
                  mvrnorm(n, mean+c(1.5*sqrt(covar[1,1]),1.0*sqrt(covar[2,2]),1.0*sqrt(covar[3,3]),1.0*sqrt(covar[4,4])), 
                          covar))
colnames(Data4)<-c("idfile","PepSeq","RT","TotalArea","MassAccu","FWHM")

#Merge all four type of disturbances + in-control observations
Data0<-rbind(Data1,Data2,Data3,Data4)
Data0<-reshape(Data0, idvar = "idfile", timevar = "PepSeq", direction = "wide")
RESPONSE<-c("NOGO")
Data0 <- cbind(Data0,RESPONSE)
Data0<-rbind(S0,Data0)

#Sample 1000 observations for training
Data<-sample_n(Data0, 10000)
#Sample 100 observations for test
Data.test<-sample_n(Data0,1000)
#Run RF model
fit <- train(as.factor(RESPONSE) ~ ., 
             data = Data[,-1], 
             method="rf",
             preProcess = c("center", "scale", "nzv"),
             tuneGrid = data.frame(mtry = 6))

#Predict RF probabilities and measure performance
Predict<-predict(fit, Data.test)
Predict.prob<-predict(fit, Data.test, type="prob")
confusionMatrix(as.factor(Data.test$RESPONSE), Predict,positive='NOGO')

#Model agnostics
importance<-varImp(fit)$importance
plot(varImp(fit))
pdpTA <- partial(fit, pred.var = "TotalArea.LVNELTEFAK",plot = TRUE, prob=TRUE)
pdpRT <- partial(fit, pred.var = "RT.LVNELTEFAK", plot = TRUE, prob=TRUE)
pdpMA <- partial(fit, pred.var = "MassAccu.LVNELTEFAK", plot = TRUE, prob=TRUE)
pdpFWHM <- partial(fit, pred.var = "FWHM.LVNELTEFAK", plot = TRUE, prob=TRUE)
grid.arrange(pdpRT, pdpTA, pdpMA, pdpFWHM, ncol = 4)
pdp1<- partial(fit, pred.var = c("TotalArea.LVNELTEFAK","RT.LVNELTEFAK"), plot = TRUE, prob=TRUE, chull=TRUE)
pdp2<- partial(fit, pred.var = c("MassAccu.LVNELTEFAK","RT.LVNELTEFAK"), plot = TRUE, prob=TRUE, chull=TRUE)
pdp3<- partial(fit, pred.var = c("FWHM.LVNELTEFAK", "RT.LVNELTEFAK"), plot = TRUE, prob=TRUE, chull=TRUE)
grid.arrange(pdp1, pdp2, pdp3,ncol = 3)
