transient_DR <- function(start_time, end_time, Final_DR, Nat_mortality, nm) {
  
  # calculate target_DR based on transient timescales
  
  # set timesteps
  years <- start_time:end_time

  # calculate moving DR array
  target_DR <- 1 - (1 - Final_DR)*(1 - exp(-1 * Nat_mortality[nm] * years))
  
  return(target_DR)
  
}