
library(h2o)
#generate in-control observations
source("sample_density_function.R")
source("add_features.R")
source("Factorial Runs/ml_algo.R")

for(sim.size in c(seq(10,1000,10)))

data.set <- data.frame()
sample_data<- list()
###### In control observation ~ 5* sim size  the of the actual j6
for(k in 1:nlevels(guide.set.scale$peptide)){ 
  sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k], sim.size*10)
  sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size*10),
                         sample_data_k,
                         RESPONSE= c("PASS"))
  sample_data[[k]] <- sample_data_k
}
data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
data <- data.set[, !duplicated(colnames(data.set))]



 for(n in c(1)){
  comb_n <-  t(combn(LETTERS[1:10],n))
  for(i in c(1:nrow(comb_n))){
    for(j in ncol(comb_n)){
      sample_data<- list()

      #small mean shift up for all peptides all metrics
      if(comb_n[i,j]=="A"){  
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                                 sample_data_k[1]+0.5*IQR(sample_data_k[,1]), 
                                 sample_data_k[2]+0.5*IQR(sample_data_k[,2]), 
                                 sample_data_k[3]+0.5*IQR(sample_data_k[,3]), 
                                 sample_data_k[4]+0.5*IQR(sample_data_k[,4]),
                                 RESPONSE= c("FAIL"))
          sample_data[[k]] <- sample_data_k
          
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large mean shift up for all peptides all metrics
      else if(comb_n[i,j]=="B"){  
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                                 sample_data_k[1]+2*IQR(sample_data_k[,1]), 
                                 sample_data_k[2]+2*IQR(sample_data_k[,2]), 
                                 sample_data_k[3]+2*IQR(sample_data_k[,3]), 
                                 sample_data_k[4]+2*IQR(sample_data_k[,4]),
                                 RESPONSE= c("FAIL"))
          sample_data[[k]] <- sample_data_k
          
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #small mean shift down for all peptides all metrics
      else if(comb_n[i,j]=="C"){  
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                                 sample_data_k[1]-0.5*IQR(sample_data_k[,1]), 
                                 sample_data_k[2]-0.5*IQR(sample_data_k[,2]), 
                                 sample_data_k[3]-0.5*IQR(sample_data_k[,3]), 
                                 sample_data_k[4]-0.5*IQR(sample_data_k[,4]),
                                 RESPONSE= c("FAIL"))
          sample_data[[k]] <- sample_data_k
          
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large mean shift down for all peptides all metrics
      else if(comb_n[i,j]=="D"){  
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          sample_data_k <- cbind(peptide=rep(levels(guide.set.scale$peptide)[k],sim.size),
                                 sample_data_k[1]-2*IQR(sample_data_k[,1]), 
                                 sample_data_k[2]-2*IQR(sample_data_k[,2]), 
                                 sample_data_k[3]-2*IQR(sample_data_k[,3]), 
                                 sample_data_k[4]-2*IQR(sample_data_k[,4]),
                                 RESPONSE= c("FAIL"))
          sample_data[[k]] <- sample_data_k
          
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample_data))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #small mean drift up for all peptides all metrics
      else if(comb_n[i,j]=="E"){  
        sample <- list()
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[k],1),
                                                           sample_data[m,1]+log(m,base=2)*IQR(sample_data[,1]), 
                                                           sample_data[m,2]+log(m,base=2)*IQR(sample_data[,2]),
                                                           sample_data[m,3]+log(m,base=2)*IQR(sample_data[,3]),
                                                           sample_data[m,4]+log(m,base=2)*IQR(sample_data[,4]),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", paste("RT",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("TotalArea",levels(guide.set.scale$peptide)[k],sep= "."),
                                       paste("MassAccu",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("FWHM",levels(guide.set.scale$peptide)[k],sep= "."),
                                       "RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large mean drift up for all peptides all metrics
      else if(comb_n[i,j]=="F"){  
        sample <- list()
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[k],1),
                                                           sample_data[m,1]+log(m,base=10)*IQR(sample_data[,1]), 
                                                           sample_data[m,2]+log(m,base=10)*IQR(sample_data[,2]),
                                                           sample_data[m,3]+log(m,base=10)*IQR(sample_data[,3]),
                                                           sample_data[m,4]+log(m,base=10)*IQR(sample_data[,4]),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", paste("RT",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("TotalArea",levels(guide.set.scale$peptide)[k],sep= "."),
                                       paste("MassAccu",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("FWHM",levels(guide.set.scale$peptide)[k],sep= "."),
                                       "RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #small mean drift down for all peptides all metrics
      else if(comb_n[i,j]=="G"){  
        sample <- list()
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[k],1),
                                                           sample_data[m,1]-log(m,base=2)*IQR(sample_data[,1]), 
                                                           sample_data[m,2]-log(m,base=2)*IQR(sample_data[,2]),
                                                           sample_data[m,3]-log(m,base=2)*IQR(sample_data[,3]),
                                                           sample_data[m,4]-log(m,base=2)*IQR(sample_data[,4]),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", paste("RT",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("TotalArea",levels(guide.set.scale$peptide)[k],sep= "."),
                                       paste("MassAccu",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("FWHM",levels(guide.set.scale$peptide)[k],sep= "."),
                                       "RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large mean drift down for all peptides all metrics
      else if(comb_n[i,j]=="H"){  
        sample <- list()
        for(k in 1:nlevels(guide.set.scale$peptide)){ 
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k],sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[k],1),
                                                           sample_data[m,1]-log(m,base=10)*IQR(sample_data[,1]), 
                                                           sample_data[m,2]-log(m,base=10)*IQR(sample_data[,2]),
                                                           sample_data[m,3]-log(m,base=10)*IQR(sample_data[,3]),
                                                           sample_data[m,4]-log(m,base=10)*IQR(sample_data[,4]),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide",paste("RT",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("TotalArea",levels(guide.set.scale$peptide)[k],sep= "."),
                                       paste("MassAccu",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("FWHM",levels(guide.set.scale$peptide)[k],sep= "."),
                                       "RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      # small Cyclic pattern
      else if(comb_n[i,j]=="I"){
        sample <- list()
        for(k in 1:nlevels(guide.set$peptide)){
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k], sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[j],1),
                                                           sample_data[m,1]+sin(2*pi*m/sim.size), 
                                                           sample_data[m,2]+sin(2*pi*m/sim.size), 
                                                           sample_data[m,3]+sin(2*pi*m/sim.size),
                                                           sample_data[m,4]+sin(2*pi*m/sim.size),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", paste("RT",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("TotalArea",levels(guide.set.scale$peptide)[k],sep= "."),
                                       paste("MassAccu",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("FWHM",levels(guide.set.scale$peptide)[k],sep= "."),
                                       "RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      #large Cyclic pattern
      else if(comb_n[i,j]=="J"){
        sample <- list()
        for(k in 1:nlevels(guide.set$peptide)){
          sample_data_k <- data.frame()
          sample_data <- sample_density(guide.set.scale,guide.set.scale$peptide[k], sim.size)
          for(m in 1:(sim.size)){
            sample_data_k <-rbind(sample_data_k,data.frame(rep(levels(guide.set.scale$peptide)[j],1),
                                                           sample_data[m,1]+2*sin(2*pi*m/sim.size), 
                                                           sample_data[m,2]+2*sin(2*pi*m/sim.size), 
                                                           sample_data[m,3]+2*sin(2*pi*m/sim.size),
                                                           sample_data[m,4]+2*sin(2*pi*m/sim.size),
                                                           "FAIL"))
          }
          colnames(sample_data_k) <- c("peptide", paste("RT",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("TotalArea",levels(guide.set.scale$peptide)[k],sep= "."),
                                       paste("MassAccu",levels(guide.set.scale$peptide)[k],sep= "."), 
                                       paste("FWHM",levels(guide.set.scale$peptide)[k],sep= "."),
                                       "RESPONSE")
          sample[[k]] <- sample_data_k
        }
        data.set <- do.call("cbind",add_features(guide.set.scale,  sample))
        data.set <- data.set[, !duplicated(colnames(data.set))]
        data <-rbind(data,data.set)
      }
      
      
    }
    
  }
  
 }



## 80% of the sample size
#smp_size <- floor(0.8 * nrow(data))

set.seed(123)
train <- data

#train_ind <- sample(seq_len(nrow(data)), size = smp_size)

#train <- data[train_ind,]
#test <- data[-train_ind,]

#launch h2o cluster
localH2O <- h2o.init(nthreads = -1)

#import r objects to h2o cloud
train_h2o <- as.h2o(train)
test_h2o <- as.h2o(test)

## run our first predictive model
rf_model <- h2o.randomForest(         ## h2o.randomForest function
  training_frame = train_h2o,        ## the H2O frame for training
  validation_frame = test_h2o,      ## the H2O frame for validation (not required)
  x= colnames(train_h2o[,c(2:5,7:162)]),
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

cf <- data.frame(h2o.confusionMatrix(rf_model))
accuracy[i] = (1- cf["Totals","Error"]) * 100 



