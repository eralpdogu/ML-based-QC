library(caret)
library(lime)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
#######SRM DATA for the paper#####################

setwd("/Users/ed/Github/ML-based-QC")
Test.set <- read.csv('test_lumos_QCloud_DDA_paper.csv')
colnames(Test.set)[2]<-"peptide"

guide.set <- read.csv('lumos_training_set.csv')

#remove repeated measurements and reshape the dataset
options(scipen = 999)
setwd("/Users/ed/Dropbox/3. MSstatsQC paper 3/4. DATA")
guide.set <-read.csv('lumos_training_set.csv')
guide.set<-guide.set[,-c(5:6)]
ind <- which(with( guide.set, (guide.set$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | guide.set$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
guide.set<-guide.set[-ind,]
guide.set$PepSeq<- gsub("\\(Carbamidomethyl\\)","",guide.set$PepSeq)
guide.set<-guide.set[,-c(2)]
colnames(guide.set)<-c('idfile', 'peptide','RT', 'TotalArea','MassAccu','FWHM')
guide.set[guide.set$RT==0,]<-NA
guide.set<-guide.set[complete.cases(guide.set),]
guide.set$peptide<-as.factor(guide.set$peptide)
guide.set$RT<-as.numeric.factor((guide.set$RT))
guide.set$TotalArea<-as.numeric(as.character(guide.set$TotalArea))
guide.set$MassAccu<-as.numeric(as.character(guide.set$MassAccu))
guide.set$FWHM<-as.numeric(as.character(guide.set$FWHM))

#remove repeated measurements and reshape the test dataset
options(scipen = 999)
setwd("/Users/ed/Dropbox/3. MSstatsQC paper 3/4. DATA")
test.set <-read.csv('lumos_all_set.csv')
ind <- which(with( test.set, (test.set$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | test.set$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
test.set<-test.set[-ind,]
test.set$PepSeq<- gsub("\\(Carbamidomethyl\\)","",test.set$PepSeq)
test.set<-test.set[,-c(2)]
colnames(test.set)<-c('idfile', 'peptide','RT', 'TotalArea', 'FWHM','MassAccu')
test.set[test.set$RT=='NULL',]<-NA
test.set<-test.set[complete.cases(test.set),]
test.set$RT<-as.numeric(as.character(test.set$RT))
test.set$TotalArea<-as.numeric(as.character(test.set$TotalArea))
test.set$MassAccu<-as.numeric(as.character(test.set$MassAccu))
test.set$FWHM<-as.numeric(as.character(test.set$FWHM))

Test.set<-test.set[1:1500,]
Test.set<-cbind(Test.set, RESPONSE=rep("NA",1500))

setwd("/Users/ed/GitHub/ML-based-QC/Factorial Runs v2")

####DDA DATA FOR PAPER#####################################
setwd("/Users/ed/Github/ML-based-QC")
guide.set <- read.csv('qtrap_training_set.csv')
colnames(guide.set)[3]<-"peptide"

#remove repeated measurements and reshape the dataset
#ind <- which(with( guide.set, (guide.set$peptide=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | guide.set$peptide=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
#guide.set<-guide.set[-ind,]
#guide.set<-guide.set[,-2]
guide.set$peptide<- gsub("\\(Carbamidomethyl\\)","",guide.set$peptide)
guide.set[guide.set$RT==0,]<-NA
guide.set[complete.cases(guide.set),]
guide.set$RT<-as.numeric.factor(guide.set$RT)
guide.set$TotalArea<-as.numeric.factor(guide.set$TotalArea)
guide.set$peptide<-as.factor(guide.set$peptide)
guide.set<-guide.set[,-c(2,4)]

guide.set.scale<-cbind(guide.set[,1:2],scale(guide.set[,-c(1:2)]))

test.set <- read.csv('qtrap_all_set.csv')
#remove repeated measurements and reshape the test dataset
options(scipen = 999)
setwd("/Users/ed/Dropbox/3. MSstatsQC paper 3/4. DATA")
test.set <-read.csv('all-lumos-2017-with-annotations.csv')
ind <- which(with( test.set, (test.set$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | test.set$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
test.set<-test.set[-ind,]
test.set$PepSeq<- gsub("\\(Carbamidomethyl\\)","",test.set$PepSeq)
test.set<-test.set[,-c(2:4, 7)]
colnames(test.set)<-c('idfile', 'peptide','RT', 'TotalArea', 'FWHM','MassAccu')
test.set[test.set$RT=='NULL',]<-NA
test.set<-test.set[complete.cases(test.set),]
test.set$RT<-as.numeric(as.character(test.set$RT))
test.set$TotalArea<-as.numeric(as.character(test.set$TotalArea))
test.set$MassAccu<-as.numeric(as.character(test.set$MassAccu))
test.set$FWHM<-as.numeric(as.character(test.set$FWHM))

#test set with annotations only
options(scipen = 999)
setwd("/Users/ed/Dropbox/3. MSstatsQC paper 3/4. DATA")
test.set <-read.csv('all-lumos-2017-with-annotations.csv')
ind <- which(with( test.set, (test.set$PepSeq=="EYEATLEEC(Carbamidomethyl)C(Carbamidomethyl)AK" | test.set$PepSeq=="TC(Carbamidomethyl)VADESHAGC(Carbamidomethyl)EK") ))
test.set<-test.set[-ind,]
test.set$PepSeq<- gsub("\\(Carbamidomethyl\\)","",test.set$PepSeq)
test.set<-test.set[test.set$Op_annotation!="NULL",]
test.set<-test.set[,-c(2:4, 7)]
colnames(test.set)<-c('idfile', 'peptide','RT', 'TotalArea', 'FWHM','MassAccu')
test.set[test.set$RT=='NULL',]<-NA
test.set<-test.set[complete.cases(test.set),]
test.set$RT<-as.numeric(as.character(test.set$RT))
test.set$TotalArea<-as.numeric(as.character(test.set$TotalArea))
test.set$MassAccu<-as.numeric(as.character(test.set$MassAccu))
test.set$FWHM<-as.numeric(as.character(test.set$FWHM))

test.set.DDA.anno<-test.set

