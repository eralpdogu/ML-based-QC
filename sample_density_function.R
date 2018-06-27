sample_density <- function(guide.set.scale,peptide, n){
  sample_data<-c()
  
  dat.dens = density(guide.set.scale[guide.set.scale$peptide == peptide,3], n=2^10)
  sim.sample.RT = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale[guide.set.scale$peptide == peptide,4], n=2^10)
  sim.sample.TotalArea = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale[guide.set.scale$peptide == peptide,5], n=2^10)
  sim.sample.MassAccu = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  dat.dens = density(guide.set.scale[guide.set.scale$peptide == peptide,6], n=2^10)
  sim.sample.FWHM = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  sample_data <- data.frame(sim.sample.RT,sim.sample.TotalArea,sim.sample.MassAccu, sim.sample.FWHM)
  names(sample_data) <- c(colnames(guide.set.scale)[3],
                          colnames(guide.set.scale)[4],
                          colnames(guide.set.scale)[5],
                          colnames(guide.set.scale)[6])
  
  return(sample_data)
}

