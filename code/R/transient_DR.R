transient_DR <- function(nat_mortality, t, time1, time2, final_DR, estimate) {
  
  # calculate target_DR based on transient timescales
  
  nm <- ifelse(estimate == 'low', 1, ifelse(estimate == 'correct', 2, 3))
  
  # set timesteps
  years <- 0:time2
  
  # set difference in M
  difference <- 0.05

  # calculate moving DR
  y <- 1 - (1 - final_DR)*(1 - exp(-1 * nat_mortality[nm] * years))
    
  # target DR depends on what time step we are on
  target_DR <- y[t - time1]
  
  return(target_DR)
  
}