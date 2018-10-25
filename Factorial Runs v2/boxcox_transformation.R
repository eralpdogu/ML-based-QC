bctrans<-function(metric){
  
  lambda.fm1 <- boxCox(metric ~ 1, family ="yjPower", plotit = FALSE)
  lambda.max <- lambda.fm1$x[which.max(lambda.fm1$y)]
  metric = yjPower(metric, lambda=lambda.max, jacobian.adjusted=FALSE)
  
}