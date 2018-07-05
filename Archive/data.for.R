#######DATA for the paper#####################

setwd("/Users/ed/Github/ML-based-QC")
Test.set <- read.csv('test_lumos_QCloud_DDA_paper.csv')
colnames(Test.set)[2]<-"peptide"

setwd("/Users/ed/Dropbox/2. MSstatsQC Paper 2-QCloud")
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
