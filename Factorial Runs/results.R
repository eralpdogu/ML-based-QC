setwd("C:/Users/aidata-1508/Documents/GitHub/ML-based-QC/Factorial Runs")
setwd("/Users/ed/GitHub/ML-based-QC/Factorial Runs")


results_NN<-read.csv("Factorial_Results_NN.txt")
results_RF<-read.csv("Factorial_Results_RF.txt")

results<-rbind(results_NN,results_RF)

g<-ggplot2::ggplot(aes(y=Accuracy,x=RUN), data=results_RF)+ geom_point()
g+xlab("Simulation Run #")
g+geom_smooth()

g<-ggplot2::ggplot(aes(y=False.positives,x=RUN), data=results_RF)+geom_point()
g+xlab("Simulation Run #")
g+xlab("False positive rate")
g+geom_smooth()

g<-ggplot2::ggplot(aes(y=False.negatives,x=RUN), data=results_RF)+geom_point()
g+xlab("Simulation Run #")
g+xlab("False negative rate")
g+geom_smooth()

df_melt <- melt(results_RF[,c(1,8:10)],id.vars ="RUN")
ggplot(df_melt, aes(RUN, variable)) + 
  geom_tile(aes(fill = value), colour = "white") +
  labs(x = "Simulation run ID",y = "Performance measures")+
  scale_fill_gradient(low = "white", high = "red",name = "") 

