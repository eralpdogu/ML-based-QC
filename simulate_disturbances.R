simulate_disturbances <- function(guide.set.scale, sim.size){

D<-rbind(simulate_drifts_partial(guide.set.scale, sim.size, 2),
                simulate_drifts_all(guide.set.scale, sim.size),
                simulate_step_shift_all(guide.set.scale, sim.size),
                simulate_step_shift_partial(guide.set.scale, sim.size, 2))

return(Data.set)
}
