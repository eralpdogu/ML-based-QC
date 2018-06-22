sample_density <- function(guide.set.scale,peptide, n){
  sample_data<-c()
  
  dat.dens = density(guide.set.scale$RT[guide.set.scale$peptide == peptide], n=2^10)
  sim.sample.RT = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale$TotalArea[guide.set.scale$peptide == peptide], n=2^10)
  sim.sample.TotalArea = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale$MassAccu[guide.set.scale$peptide == peptide], n=2^10)
  sim.sample.MassAccu = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale$FWHM[guide.set.scale$peptide == peptide], n=2^10)
  sim.sample.FWHM = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  sample_data <- data.frame(sim.sample.RT,sim.sample.TotalArea,sim.sample.MassAccu, sim.sample.FWHM)
  names(sample_data) <- c("RT","TotalArea ","MassAccu","FWHM")
  
  return(sample_data)
}

