transient_DR <- function(M, t, time1, time2, final_DR, estimate) {
  
  # calculate target_DR based on transient timescales
  
  # set timesteps
  years <- 0:time2
  
  # set difference in M
  difference <- 0.05
  
  # low M
  if (estimate == 'low') {
    m <- M - difference
    y <- 1 - (1 - final_DR)*(1 - exp(-1 * m * years))
  }
  
  # correct M
  if (estimate == 'correct') {
    y <- 1 - (1 - final_DR)*(1 - exp(-1 * M * years))
  }
  
  # high M
  if (estimate == 'high') {
    m <- M + difference
    y <- 1 - (1 - final_DR)*(1 - exp(-1 * m * years))
  }
  
  target_DR <- y[t - time1]
  
  return(target_DR)
  
}