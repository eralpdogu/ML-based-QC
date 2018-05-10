# Similarity matrix in r

#train data 
#S0 
#N0

#test data
#Data
#NW

#makes sure the data is numeric
sw <- sapply(Data[,c(-1,-50)],as.numeric)
s0 <- sapply(S0[,c(-1,-50)],as.numeric)

#install.packages("resemble")
library(resemble) # cosine
#?fDiss

#take the average representative point of train data 
avg_s0 <- colMeans(s0)

#finding cosine similarity between test data and train data points
fDiss(s0,t(as.data.frame(sw[1,])), method = "cosine",center = F, scaled = F)

fDiss(s0,t(as.data.frame(sw[1,])), method = "cosine",center = F, scaled = F)


#takes a sample of test data on find simlarity on real time basis
for(i in 1:length(Data)){
  mycosine(Data[i,],S0)
  
}
