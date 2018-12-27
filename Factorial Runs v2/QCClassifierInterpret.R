QCClassifierInterpret<-function(trainData, Test.set.scale, model, caseNumber){
  explainer <- lime(train, rf_model, n_bins = 5)
  explanation <- explain(Test.set.scale, explainer, n_bins = 5)
  explanation_caret <- explain(
    x = Test.set, 
    explainer = explainer, 
    n_permutations = 5000,
    dist_fun = "gower",
    kernel_width = .75,
    n_features = 10, 
    feature_select = "highest_weights",
    labels = "Yes"
  )
  
  plot_features(explanation[40:50,])
}