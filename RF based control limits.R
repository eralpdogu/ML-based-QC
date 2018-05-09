########Control limits#############################
N0<-length(S0[,1])
Nw<-10
p0<-c()
p1<-c()

for (i in Nw:length(S0[,1])){
  SW<-S0[(i-Nw+1):i,]
  SW$RESPONSE<-'NOGO'
  Data.model<-do.call(cbind.data.frame, Map('c',S0, SW))
  fit <- train(as.factor(RESPONSE) ~ ., 
               data = Data.model, 
               method="rf",
               preProcess = c("center", "scale", "nzv"),
               trainControl=trainControl( method="oob" ),  
               tuneGrid = data.frame(mtry = 6))
  
  Predict<-predict(fit, type='prob')
  nTrees[i]<-fit$finalModel$ntree
  p1[i]<-sum(Predict[(N0+1):(N0+Nw),2])/Nw

}

n = length(p0) 
B = 1000
result = rep(NA, B)
for (i in 1:B) {
  boot.sample = sample(n, replace = TRUE)
  result[i] = quantile(p1[boot.sample],0.9973,na.rm = T)
}
CL<-mean(result)
