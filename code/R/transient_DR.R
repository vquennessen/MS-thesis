transient_DR <- function(nat_mortality, t, time1, time2, final_DR, nm) {
  
  # calculate target_DR based on transient timescales
  
  # set timesteps
  years <- 0:time2

  # calculate moving DR
  y <- 1 - (1 - final_DR)*(1 - exp(-1 * nat_mortality[nm] * years))
    
  # target DR depends on what time step we are on
  target_DR <- y[t - time1]
  
  return(target_DR)
  
}