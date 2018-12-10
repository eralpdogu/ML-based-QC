
# Step shift --------------------------------------------------------------

QcClassifier_data_step <- function(data,nmetrics,factor.names,sim.size,peptide.colname){
  
  if(!is.factor(new_data[,paste(peptide.colname)])){
    new_data$peptide <-as.factor(new_data$peptide.colname)  
  }
  #factorial matrix  
  factorial <- FrF2(2^nmetric, nmetric,factor.names=factor.names)
  
  
  tag_neg <- 0
  data <- data.frame(NULL)
  
  for(i in 1:nrow(factorial)){
    data.set <- data.frame(NULL)
    if(all(factorial[i,]== rep(-1,nmetric))){
      ####### In cntrol observation ~ 5* sim size  the of the actual 
      sample_data_k <- sample_density(new_data, sim.size*(2^(nmetric)-1))
    }
    
    else{
      ###### Base Data set to begin with 
      sample_data_k <- sample_density(new_data, sim.size)
      #sample_data_k <- robust.scale(sample_data_k)
      
      for(j in 1:ncol(sample_data_k)){
        #change in a metric for some peptides
        if(factorial[i,j]== "1" & colnames(factorial[i,j])==colnames(sample_data_k)[j]){ 
          beta=runif(sim.size,-5,5)
          sample_data_k[,j] <- sample_data_k[,j] + beta*mad(sample_data_k[,j])
          tag_neg <- 1 
          
        }# column ends 
      }
    }
    data.set <- rbind(data.set,add_features(sample_data_k))
    #data.set[,"peptide"] <- NULL 
    if(tag_neg == 1){
      data.set$RESPONSE <- c("FAIL")
      tag_neg <- 0
    }
    else{
      data.set$RESPONSE <- c("PASS")
    }
    data <- data[,order(names(data))]
    data.set <- data.set[,order(names(data.set))]
    data <-rbind(data,data.set)
  }
  
  data <- data[sample(nrow(data), nrow(data)), ] # shuffle the data
  data$RESPONSE <- as.factor(data$RESPONSE)
  
  return(data) 
}

# Variance change ---------------------------------------------------------

QcClassifier_data_var <- function(data,nmetrics,factor.names,sim.size,peptide.colname){
  
  if(!is.factor(new_data[,paste(peptide.colname)])){
    new_data$peptide <-as.factor(new_data$peptide.colname)  
  }
  #factorial matrix  
  factorial <- FrF2(2^nmetric, nmetric,factor.names=factor.names)
  
  
  tag_neg <- 0
  data <- data.frame(NULL)
  
  for(i in 1:nrow(factorial)){
    data.set <- data.frame(NULL)
    if(all(factorial[i,]== rep(-1,nmetric))){
      ####### In cntrol observation ~ 5* sim size  the of the actual 
      sample_data_k <- sample_density(new_data, sim.size*(2^(nmetric)-1))
    }
    
    else{
      ###### Base Data set to begin with 
      sample_data_k <- sample_density(new_data, sim.size)
      #sample_data_k <- robust.scale(sample_data_k)
      
      for(j in 1:ncol(sample_data_k)){
        #change in a metric for some peptides
        if(factorial[i,j]== "1" & colnames(factorial[i,j])==colnames(sample_data_k)[j]){ 
          beta=runif(sim.size,-5,5)
          sample_data_k[,j] <- sample_data_k[,j]*beta
          tag_neg <- 1 
          
        }# column ends 
      }
    }
    data.set <- rbind(data.set,add_features(sample_data_k))
    #data.set[,"peptide"] <- NULL 
    if(tag_neg == 1){
      data.set$RESPONSE <- c("FAIL")
      tag_neg <- 0
    }
    else{
      data.set$RESPONSE <- c("PASS")
    }
    data <- data[,order(names(data))]
    data.set <- data.set[,order(names(data.set))]
    data <-rbind(data,data.set)
  }
  
  data <- data[sample(nrow(data), nrow(data)), ] # shuffle the data
  data$RESPONSE <- as.factor(data$RESPONSE)
  
  return(data) 
}


# Linear drift ------------------------------------------------------------

QcClassifier_data_linear <- function(data,nmetrics,factor.names,sim.size,peptide.colname){
  
  if(!is.factor(new_data[,paste(peptide.colname)])){
    new_data$peptide <-as.factor(new_data$peptide.colname)  
  }
  #factorial matrix  
  factorial <- FrF2(2^nmetric, nmetric,factor.names=factor.names)
  
  
  tag_neg <- 0
  data <- data.frame(NULL)
  
  for(i in 1:nrow(factorial)){
    data.set <- data.frame(NULL)
    if(all(factorial[i,]== rep(-1,nmetric))){
      ####### In cntrol observation ~ 5* sim size  the of the actual 
      sample_data_k <- sample_density(new_data, sim.size*(2^(nmetric)-1))
    }
    
    else{
      ###### Base Data set to begin with 
      sample_data_k <- sample_density(new_data, sim.size)
      #sample_data_k <- robust.scale(sample_data_k)
      
      for(j in 1:ncol(sample_data_k)){
        #change in a metric for some peptides
        if(factorial[i,j]== "1" & colnames(factorial[i,j])==colnames(sample_data_k)[j]){ 
          beta=runif(sim.size,0,2)
          for(i in 1:sim.size){
          beta[i]=beta[i]*(i-sim.size)}
          sample_data_k[,j] <- sample_data_k[,j] + beta*mad(sample_data_k[,j])
          tag_neg <- 1 
          
        }# column ends 
      }
    }
    data.set <- rbind(data.set,add_features(sample_data_k))
    #data.set[,"peptide"] <- NULL 
    if(tag_neg == 1){
      data.set$RESPONSE <- c("FAIL")
      tag_neg <- 0
    }
    else{
      data.set$RESPONSE <- c("PASS")
    }
    data <- data[,order(names(data))]
    data.set <- data.set[,order(names(data.set))]
    data <-rbind(data,data.set)
  }
  
  data <- data[sample(nrow(data), nrow(data)), ] # shuffle the data
  data$RESPONSE <- as.factor(data$RESPONSE)
  
  return(data) 
}
