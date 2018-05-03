########Control limits#############################
N0<-length(S0[,1])
Nw<-10
p0<-c()
p1<-c()

for (i in N0+1:length(S0[,1])){
  SW<-S0
  SW$RESPONSE<-'NOGO'
  Data.model<-do.call(cbind.data.frame, Map('c',S0, SW))
  
  fit <- train(as.factor(RESPONSE) ~ ., data = Data.model, method="rf",  preProcess = c("center", "scale", "nzv")
               , trainControl=trainControl( method="oob" ))
  
  Predict<-predict(fit, type='prob')
  p0[i]<-sum(Predict[,1])/length(S0[,1])

}

n = length(p0) 
B = 1000
result = rep(NA, B)
for (i in 1:B) {
  boot.sample = sample(n, replace = TRUE)
  result[i] = quantile(p0[boot.sample],0.997,na.rm = T)
}
CL<-mean(result)
