sample_density <- function(guide.set, n){
  sample_data<-c()
  guide.set.scale<-data.frame(NULL)
  for(i in 1:nlevels(guide.set$peptide)){
  guide.set.temp<-robust.scale(guide.set[guide.set$peptide==levels(guide.set$peptide)[i],3:6])
  guide.set.scale<-rbind(guide.set.scale, guide.set.temp)
  }
  guide.set.scale<-cbind(guide.set[,1:2],guide.set.scale)
  #dat.dens = density(guide.set[guide.set$peptide == peptide,3], n=2^10)
  dat.dens = density(guide.set.scale[1:8,3], n=2^10)
  sim.sample.RT = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,4], n=2^10)
  dat.dens = density(guide.set.scale[,4], n=2^10)
  sim.sample.TotalArea = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,5], n=2^10)
  dat.dens = density(guide.set.scale[,5], n=2^10)
  sim.sample.MassAccu = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,6], n=2^10)
  dat.dens = density(guide.set.scale[,6], n=2^10)
  sim.sample.FWHM = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  sample_data <- data.frame(sim.sample.RT,sim.sample.TotalArea,sim.sample.MassAccu, sim.sample.FWHM)
  names(sample_data) <- c("RT", "TotalArea", "MassAccu", "FWHM")
  
  return(sample_data)
}

