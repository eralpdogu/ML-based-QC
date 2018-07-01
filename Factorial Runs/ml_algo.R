
ml_algo <- function(data,num_run){
  
  ## 80% of the sample size
  smp_size <- floor(0.8 * nrow(data))
  
  set.seed(123)
  train_ind <- sample(seq_len(nrow(data)), size = smp_size)
  
  train <- data[train_ind,]
  test <- data[-train_ind,]
  
  
  dir.create(paste("Factorial_Rf_Runs/",num_run,sep = ""))
  
  
  
  # #Random Forests ###########################################
  # print("Random Forest: Train Data using all features :")
  # fit_all <- train(y=train[,"RESPONSE"],x=subset(train,select = -c(RESPONSE,peptide)),
  #                  method="rf",
  #                  preProcess = c("center", "scale", "nzv"),
  #                  tuneGrid = data.frame( mtry=floor(sqrt(ncol(train))) )) # change mtry as square root of number of predictors
  # 
  # print(fit_all$results)
  # 
  # Predict<-predict(fit_all, test)
  # print("Random Forest: Test Data using all features :")
  # cf <- confusionMatrix(as.factor(train$RESPONSE), Predict,positive='FAIL')
  # 
  # save(fit_all,cf,file = paste("Factorial_Rf_Runs/",num_run,"/","Random_Forest_Results",num_run,sep = ""))
  # 
  # return(fit_all)
  
  
  #launch h2o cluster
  localH2O <- h2o.init(nthreads = -1)

  #import r objects to h2o cloud
  train_h2o <- as.h2o(train)
  test_h2o <- as.h2o(test)
  
  ## run our first predictive model
  rf_model <- h2o.randomForest(         ## h2o.randomForest function
    training_frame = train_h2o,        ## the H2O frame for training
    validation_frame = test_h2o,      ## the H2O frame for validation (not required)
    x= colnames(train_h2o[,c(2:5,7:102)]),
    y= "RESPONSE",
    model_id = "rf_model",    ## name the model in H2O
    ntrees = 200,                  ##   not required, but helps use Flow
                                   ## use a maximum of 200 trees to create the
                                   ##  random forest model. The default is 50.
    stopping_rounds = 2,           ## Stop fitting new trees when the 2-tree
                                   ##  average is within 0.001 (default) of 
                                   ##  the prior two 2-tree averages.
                                   ##  Can be thought of as a convergence setting
    score_each_iteration = T,     
    seed = 123) 
  
  
  ###############################################################################
  summary(rf_model)  
  
  h2o.saveModel(rf_model, path = paste("Factorial_Rf_Runs/",num_run,sep = ""))
  png(filename=paste("Factorial_Rf_Runs/",num_run,"/",num_run,".png",sep = ""))
  plot(rf_model)

  dev.off()
  
  h2o.shutdown(prompt = F)
  
  return(rf_model)
  
  
  # #Neural Network ################################################################################
  # #launch h2o cluster
  # localH2O <- h2o.init(nthreads = -1)
  # 
  # #import r objects to h2o cloud
  # train_h2o <- as.h2o(train)
  # test_h2o <- as.h2o(test)
  # 
  # dl_model <- h2o.deeplearning(
  #   model_id="dl_model", 
  #   training_frame=train_h2o, 
  #   validation_frame = test_h2o,
  #   x= colnames(train_h2o[,c(2:5,7:102)]),
  #   y= "RESPONSE",
  #   activation ="Tanh",  
  #   hidden=c(20,20), 
  #   standardize = TRUE, #standardizes the data
  #   loss= "CrossEntropy",
  #   stopping_metric="logloss",
  #   stopping_rounds = 10,
  #   stopping_tolerance=0.001,
  #   adaptive_rate = TRUE,
  #   shuffle_training_data = TRUE, 
  #   rate = 0.005, # Defaults to 0.005 adaptive enabled so cannot specify the learning rare 
  #   mini_batch_size = 1,# defaults to 1 
  #   epochs=200,
  #   seed = 123, # give seed 
  #   export_weights_and_biases = T, # export weights and biases defaults to false
  #   reproducible = T # Force reproducibility on small data (will be slow - only uses 1 thread). Defaults to FALSE.
  # )
  # print(summary(dl_model))
  # 
  # 
  # 
  # h2o.saveModel(dl_model, path = paste("Factorial Runs/",num_run,sep = ""))
  # png(filename=paste("Factorial Runs/",num_run,"/",num_run,".png",sep = ""))
  # plot(dl_model)
  # 
  # dev.off()
  # 
  # 
  # 
  # h2o.shutdown(prompt = F)
  # 
  # return(dl_model)
  
}

