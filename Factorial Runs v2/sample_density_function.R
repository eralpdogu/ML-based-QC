sample_density <- function(guide.set, n){
  sample_data<-c()
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,3], n=2^10)
  dat.dens = density(guide.set[1:8,3], n=2^10)
  sim.sample.RT = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,4], n=2^10)
  dat.dens = density(guide.set[,4], n=2^10)
  sim.sample.TotalArea = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,5], n=2^10)
  dat.dens = density(guide.set[,5], n=2^10)
  sim.sample.MassAccu = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,6], n=2^10)
  dat.dens = density(guide.set[,6], n=2^10)
  sim.sample.FWHM = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  sample_data <- data.frame(sim.sample.RT,sim.sample.TotalArea,sim.sample.MassAccu, sim.sample.FWHM)
  names(sample_data) <- c("RT", "TotalArea", "MassAccu", "FWHM")
  
  return(sample_data)
}

