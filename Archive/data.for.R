library(caret)
library(lime)

#######DATA for the paper#####################

setwd("/Users/ed/Github/ML-based-QC")
Test.set <- read.csv('test_lumos_QCloud_DDA_paper.csv')
colnames(Test.set)[2]<-"peptide"

setwd("C:/Users/aidata-1508/Dropbox/2. MSstatsQC Paper 2-QCloud")

guide.set <- read.csv('lumos_training_set.csv')

#remove repeated measurements and reshape the dataset
ind <- which(with( Train, (guide.set$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | guide.set$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
guide.set<-guide.set[-ind,]
guide.set<-guide.set[,-2]
guide.set$PepSeq<- gsub("\\(Carbamidomethyl\\)","",guide.set$PepSeq)
guide.set[guide.set$RT==0,]<-NA
guide.set[complete.cases(guide.set),]

guide.set.scale<-cbind(guide.set[,1:2],scale(guide.set[,-c(1:2)]))

#remove repeated measurements and reshape the test dataset
options(scipen = 999)
setwd("/Users/ed/Dropbox/3. MSstatsQC paper 3/4. DATA")
test.set <-read.csv('lumos_all_set.csv')
ind <- which(with( test.set, (test.set$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | test.set$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
test.set<-test.set[-ind,]
test.set<-test.set[,-c(5,6)]
test.set<-test.set[,-2]
test.set$PepSeq<- gsub("\\(Carbamidomethyl\\)","",test.set$PepSeq)
test.set[test.set$rt.sec.==0,]<-NA
test.set<-test.set[complete.cases(test.set),]
colnames(test.set)<-c('idfile', 'peptide','RT', 'TotalArea', 'MassAccu','FWHM')
Test.set<-test.set

setwd("/Users/ed/GitHub/ML-based-QC/Factorial Runs v2")

