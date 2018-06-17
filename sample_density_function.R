sample_density <- function(guide.set.scale,peptide, n){
  
  dat.dens = density(guide.set.scale$RT[guide.set.scale$Precursor == peptide], n=2^10)
  sim.sample.RT = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale$TotalArea[guide.set.scale$Precursor == peptide], n=2^10)
  sim.sample.TotalArea = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale$MassAccu[guide.set.scale$Precursor == peptide], n=2^10)
  sim.sample.MassAccu = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale$FWHM[guide.set.scale$Precursor == peptide], n=2^10)
  sim.sample.FWHM = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  sample_data <- list(sim.sample.RT,sim.sample.TotalArea,sim.sample.MassAccu, sim.sample.FWHM)
  names(sample_data) <- c("sim.sample.RT","sim.sample.TotalArea ","sim.sample.MassAccu","sim.sample.FWHM")
  
  return(sample_data)
}

