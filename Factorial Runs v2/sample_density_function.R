sample_density <- function(guide.set, n){
  sample_data<-c()
  guide.set.scaled<-data.frame(NULL)
  for(i in 1:nlevels(guide.set$peptide)){
  guide.set.temp<-robust.scale(guide.set[guide.set$peptide==levels(guide.set$peptide)[i],3:6])
  guide.set.scaled<-rbind(guide.set.scaled, guide.set.temp)
  }
  guide.set.scaled<-cbind(guide.set[,1:2],bctrans(guide.set.scaled$RT),
                                         bctrans(guide.set.scaled$TotalArea),
                                         bctrans(guide.set.scaled$MassAccu),
                                         bctrans(guide.set.scaled$FWHM))
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,3], n=2^10)
  dat.dens = density(guide.set.scaled[,3], n=2^10)
  sim.sample.RT = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,4], n=2^10)
  dat.dens = density(guide.set.scaled[,4], n=2^10)
  sim.sample.TotalArea = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,5], n=2^10)
  dat.dens = density(guide.set.scaled[,5], n=2^10)
  sim.sample.MassAccu = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  #dat.dens = density(guide.set[guide.set$peptide == peptide,6], n=2^10)
  dat.dens = density(guide.set.scaled[,6], n=2^10)
  sim.sample.FWHM = sample(dat.dens$x, n, replace=TRUE, prob=dat.dens$y)
  
  sample_data <- data.frame(sim.sample.RT,sim.sample.TotalArea,sim.sample.MassAccu, sim.sample.FWHM)
  names(sample_data) <- c("RT", "TotalArea", "MassAccu", "FWHM")
  
  return(sample_data)
}

