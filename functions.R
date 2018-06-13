library(caret)
library(h2o)
library(MASS)


#n = 1 for RF , 2 for ANN 
modeling <- function(n){
  if(n==1){
    #load the test and train data 
    load("Sim_Test_Train_Data.RData")
    
    #RF model
    fit <- train(RESPONSE ~ ., 
                 data = train_sim_all[,-c(1,2)],  
                 method="rf",
                 preProcess = c("center", "scale", "nzv"),
                 tuneGrid = data.frame(mtry = 6)) # change mtry as square root of number of predictors
    print(summary(fit))
    
    #Predict RF probabilities and measure performance
    Predict<-predict(fit, test_sim_all[,-c(1,2,7)]) # change the features here 
    Predict.prob<-predict(fit, test_sim_all, type="prob")
    print("Test Data",confusionMatrix(test_sim_all$RESPONSE, Predict))
    
    #Model agnostics
    importance<-varImp(fit)$importance
    plot(varImp(fit))
    
    return(fit) # return the model object 
  }
  
  else if(n==2){
    #ANN Model
    #Building Neural Network
    library(h2o)
    #generate same set of random numbers (for reproducibility)
    set.seed(121)
    
    #launch h2o cluster
    localH2O <- h2o.init(nthreads = -1)
    
    
    #import r objects to h2o cloud
    train_h2o <- as.h2o(train_sim_all)
    test_h2o <- as.h2o(test_sim_all)
    
    set.seed(100)
    
    dl_model_sim1 <- h2o.deeplearning(
      model_id="dl_model_first", 
      training_frame=train_h2o, 
      validation_frame = test_h2o,
      x= colnames(train_h2o[,3:6]),
      y= "RESPONSE",
      activatio="Tanh",  
      hidden=c(20,20), 
      standardize = TRUE, #standardizes the data
      loss= "CrossEntropy",
      stopping_metric="logloss",
      stopping_rounds = 10,
      stopping_tolerance=0.00001,
      adaptive_rate = TRUE,
      shuffle_training_data = TRUE, 
      rate = 0.005, # Defaults to 0.005 adaptive enabled so cannot specify the learning rare 
      mini_batch_size = 1,# defaults to 1 
      epochs=200,
      seed = 123, # give seed 
      export_weights_and_biases = T, # export weights and biases defaults to false
      reproducible = T # Force reproducibility on small data (will be slow - only uses 1 thread). Defaults to FALSE.
    )
    summary(dl_model_sim1)
    plot(dl_model_sim1)
    h2o.shutdown()
    
    return(dl_model_sim1) # return the model object 
    
  }
    else{
      
      Print("Illegal Response")
    }

}

# calling the function
# 1 for RF 
modeling(1) 

# calling the function
# 2 for Neural Nwtworks 
modeling(2)

