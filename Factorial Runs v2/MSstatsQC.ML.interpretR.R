QCClassifierInterpret<-function(trainData, Test.set.scale, model, caseNumber){
  
  explainer <- lime(train[,-11], rf_model, n_bins = 5)
  
  explanation_caret <- explain(
    x = Test.set.scale, 
    explainer = explainer, 
    n_labels = 1,
    n_features = 3
  )
  
  plot_features(explanation_caret[explanation_caret$case==236:240,])
}