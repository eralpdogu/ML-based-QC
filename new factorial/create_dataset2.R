library(readxl)
library(h2o)
library(caret)
library(MASS)
library(ggplot2)

factorial <- read_xlsx("factorial.xlsx",sheet = 1)

source("add_features.R")
source("ml_algo.R")
data.set <- data.frame()


acc <- data.frame(matrix(nrow=4,ncol=10))
row.names(acc) <- c(10, 25, 50, 500)
names(acc)<- seq(10,100,10)


fp <- data.frame(matrix(nrow=4,ncol=10))
row.names(fp) <- c(10, 25, 50, 500)
names(fp)<- seq(10,100,10)


fn <- data.frame(matrix(nrow=4,ncol=10))
row.names(fn) <- c(10, 25, 50, 500)
names(fn)<- seq(10,100,10)

for (sim.size in c(10, 25, 50, 500)){
  for(guide_percent_row in seq(10,100,10)){
    num_peptides <- nlevels(guide.set.scale$peptide)
    for(i in c(1)){
      sample_data<- list()
      guide.ss <- guide.set.scale[1:floor(guide_percent_row*nrow(guide.set.scale)/100),]
      ###### In control observation ~ 5* sim size  the of the actual 
      for(k in 1:nlevels(guide.set.scale$peptide)){ 
        sample_data_k <- sample_density( guide.ss, guide.ss$peptide[k], sim.size*5)
        sample_data_k <- cbind(peptide=rep(levels( guide.ss$peptide)[k],sim.size*5),
                               sample_data_k,
                               RESPONSE= c("PASS"))
        sample_data[[k]] <- sample_data_k
      }
      data.set <- do.call("cbind",add_features( guide.ss,  sample_data))
      data <- data.set[, !duplicated(colnames(data.set))]
      
      
      for(j in 2:6){
        #mean shift up for all peptides all metrics
        if(colnames(factorial[i,j])=="A"){  
          for(k in 1:nlevels( guide.ss$peptide)){ 
            sample_data_k <- sample_density( guide.ss, guide.ss$peptide[k],sim.size)
            sample_data_k <- cbind(peptide=rep(levels( guide.ss$peptide)[k],sim.size),
                                   sample_data_k[1]+runif(1,0.5,2)*IQR(sample_data_k[,1]), 
                                   sample_data_k[2]+runif(1,0.5,2)*IQR(sample_data_k[,2]), 
                                   sample_data_k[3]+runif(1,0.5,2)*IQR(sample_data_k[,3]), 
                                   sample_data_k[4]+runif(1,0.5,2)*IQR(sample_data_k[,4]),
                                   RESPONSE= c("FAIL"))
            sample_data[[k]] <- sample_data_k
            
          }
          data.set <- do.call("cbind",add_features( guide.ss,  sample_data))
          data.set <- data.set[, !duplicated(colnames(data.set))]
          data <-rbind(data,data.set)
        }
        
        
        #mean shift down for all peptides all metrics
        else if(colnames(factorial[i,j])=="B"){  
          for(k in 1:nlevels( guide.ss$peptide)){ 
            sample_data_k <- sample_density( guide.ss, guide.ss$peptide[k],sim.size)
            sample_data_k <- cbind(peptide=rep(levels( guide.ss$peptide)[k],sim.size),
                                   sample_data_k[1]-runif(1,0.5,2)*IQR(sample_data_k[,1]), 
                                   sample_data_k[2]-runif(1,0.5,2)*IQR(sample_data_k[,2]), 
                                   sample_data_k[3]-runif(1,0.5,2)*IQR(sample_data_k[,3]), 
                                   sample_data_k[4]-runif(1,0.5,2)*IQR(sample_data_k[,4]),
                                   RESPONSE= c("FAIL"))
            sample_data[[k]] <- sample_data_k
            
          }
          data.set <- do.call("cbind",add_features( guide.ss,  sample_data))
          data.set <- data.set[, !duplicated(colnames(data.set))]
          data <-rbind(data,data.set)
        }
        
        
        #mean drift up for all peptides all metrics
        else if(colnames(factorial[i,j])=="C"){  
          sample <- list()
          for(k in 1:nlevels( guide.ss$peptide)){ 
            sample_data_k <- data.frame()
            sample_data <- sample_density( guide.ss, guide.ss$peptide[k],sim.size)
            for(m in 1:(sim.size)){
              sample_data_k <-rbind(sample_data_k,data.frame(rep(levels( guide.ss$peptide)[k],1),
                                                             sample_data[m,1]+log(m,base=runif(1,2,10))*IQR(sample_data[,1]), 
                                                             sample_data[m,2]+log(m,base=runif(1,2,10))*IQR(sample_data[,2]),
                                                             sample_data[m,3]+log(m,base=runif(1,2,10))*IQR(sample_data[,3]),
                                                             sample_data[m,4]+log(m,base=runif(1,2,10))*IQR(sample_data[,4]),
                                                             "FAIL"))
            }
            colnames(sample_data_k) <- c("peptide", paste("RT",levels( guide.ss$peptide)[k],sep= "."), 
                                         paste("TotalArea",levels( guide.ss$peptide)[k],sep= "."),
                                         paste("MassAccu",levels( guide.ss$peptide)[k],sep= "."), 
                                         paste("FWHM",levels( guide.ss$peptide)[k],sep= "."),
                                         "RESPONSE")
            sample[[k]] <- sample_data_k
          }
          data.set <- do.call("cbind",add_features( guide.ss,  sample))
          data.set <- data.set[, !duplicated(colnames(data.set))]
          data <-rbind(data,data.set)
        }
        
        
        #mean drift down for all peptides all metrics
        else if(colnames(factorial[i,j])=="D"){  
          sample <- list()
          for(k in 1:nlevels( guide.ss$peptide)){ 
            sample_data_k <- data.frame()
            sample_data <- sample_density( guide.ss, guide.ss$peptide[k],sim.size)
            for(m in 1:(sim.size)){
              sample_data_k <-rbind(sample_data_k,data.frame(rep(levels( guide.ss$peptide)[k],1),
                                                             sample_data[m,1]-log(m,base=runif(1,2,10))*IQR(sample_data[,1]), 
                                                             sample_data[m,2]-log(m,base=runif(1,2,10))*IQR(sample_data[,2]),
                                                             sample_data[m,3]-log(m,base=runif(1,2,10))*IQR(sample_data[,3]),
                                                             sample_data[m,4]-log(m,base=runif(1,2,10))*IQR(sample_data[,4]),
                                                             "FAIL"))
            }
            colnames(sample_data_k) <- c("peptide", paste("RT",levels( guide.ss$peptide)[k],sep= "."), 
                                         paste("TotalArea",levels( guide.ss$peptide)[k],sep= "."),
                                         paste("MassAccu",levels( guide.ss$peptide)[k],sep= "."), 
                                         paste("FWHM",levels( guide.ss$peptide)[k],sep= "."),
                                         "RESPONSE")
            sample[[k]] <- sample_data_k
          }
          data.set <- do.call("cbind",add_features( guide.ss,  sample))
          data.set <- data.set[, !duplicated(colnames(data.set))]
          data <-rbind(data,data.set)
        }
        
        
        #Cyclic pattern
        else if(colnames(factorial[i,j])=="E"){
          sample <- list()
          for(k in 1:nlevels( guide.ss$peptide)){
            sample_data_k <- data.frame()
            sample_data <- sample_density( guide.ss, guide.ss$peptide[k], sim.size)
            for(m in 1:(sim.size)){
              sample_data_k <-rbind(sample_data_k,data.frame(rep(levels( guide.ss$peptide)[j],1),
                                                             sample_data[m,1]+runif(1,1,2)*sin(2*pi*m/sim.size), 
                                                             sample_data[m,2]+runif(1,1,2)*sin(2*pi*m/sim.size), 
                                                             sample_data[m,3]+runif(1,1,2)*sin(2*pi*m/sim.size),
                                                             sample_data[m,4]+runif(1,1,2)*sin(2*pi*m/sim.size),
                                                             "FAIL"))
            }
            colnames(sample_data_k) <- c("peptide", paste("RT",levels( guide.ss$peptide)[k],sep= "."), 
                                         paste("TotalArea",levels( guide.ss$peptide)[k],sep= "."),
                                         paste("MassAccu",levels( guide.ss$peptide)[k],sep= "."), 
                                         paste("FWHM",levels( guide.ss$peptide)[k],sep= "."),
                                         "RESPONSE")
            sample[[k]] <- sample_data_k
          }
          data.set <- do.call("cbind",add_features( guide.ss,  sample))
          data.set <- data.set[, !duplicated(colnames(data.set))]
          data <-rbind(data,data.set)
        }
        
        
      }# column ends 
      data <- data[sample(nrow(data), nrow(data)), ] # shuffle the data
      
      
      
      dl_model <- ml_algo(data,i)
      
      
      
      cf<- data.frame(h2o.confusionMatrix(dl_model,valid = T),stringsAsFactors = F)
      acc[paste(sim.size),paste(guide_percent_row)] <- (1-cf[3,3])*100
      fp[paste(sim.size),paste(guide_percent_row)] <- cf[2,3]
      fn[paste(sim.size),paste(guide_percent_row)] <- cf[1,3]
      #write.csv(factorial, file = paste("factorial_guideset_",nrow(guide.set.scale),"_",sim.size,".csv",sep = ""))
      
    }
    
  }
}

acc$sim_size <- row.names(acc)
df_melt <- melt(acc,id.vars = "sim_size")
ggplot(df_melt, aes(sim_size, variable)) + 
  geom_tile(aes(fill = value), colour = "white") +
  labs(x = "Simulation Size(10x)",y = "Percentage of rows in Guide Set")+
  scale_fill_gradient(low = "white", high = "blue",name = "Accuracy") 

fp$sim_size <- row.names(fp)
df_melt <- melt(fp,id.vars = "sim_size")
ggplot(df_melt, aes(sim_size, variable)) + 
  geom_tile(aes(fill = value), colour = "white") +
  labs(x = "Simulation Size(10x)",y = "Percentage of rows in Guide Set")+
  scale_fill_gradient(low = "white", high = "blue",name = "False Positives") 

fn$sim_size <- row.names(fn)
df_melt <- melt(fp,id.vars = "sim_size")
ggplot(df_melt, aes(sim_size, variable)) + 
  geom_tile(aes(fill = value), colour = "white") +
  labs(x = "Simulation Size(10x)",y = "Percentage of rows in Guide Set")+
  scale_fill_gradient(low = "white", high = "blue",name = "False Negatives") 


