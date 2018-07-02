setwd("C:/Users/aidata-1508/Documents/GitHub/ML-based-QC/Factorial Runs")

results_NN<-read.csv("Factorial_Results_RF.txt")
results_RF<-read.csv("Factorial_Results_NN.txt")

results<-rbind(results_NN,results_RF)

ggplot2::ggplot(aes(y=Overall.error.rate,x=Run.label, col=Method), data=results)+geom_point()
ggplot2::ggplot(aes(y=False.positives,x=Run.label, col=Method), data=results)+geom_point()
ggplot2::ggplot(aes(y=False.negatives,x=Run.label, col=Method), data=results)+geom_point()

