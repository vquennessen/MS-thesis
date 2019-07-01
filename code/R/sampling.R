sampling <- function(a, t, y, r, D, abundance_all, abundance_mature, transects, 
                     x, count_sp, nuS) {
  
  # Calculate delta
  # Based on Babcock & MacCall (2011): Eq. (13)
  delta <- r / D
  
  # Calculate probability of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  p_all <-  delta * abundance_all[a, t] * A / R0
  p_mature <-  delta * abundance_mature[a, t] * A / R0

  # TODO: figure out how to properly calculate probability
  
  # Determine if species is seen at least once
  # Dimensions = 1 * transects
  presence_all <- rbinom(transects, 1, p_all)
  presence_mature <- rbinom(transects, 1, p_mature)
  
  # Based on Babcock & MacCall (2011): Eq. (16)
  gamma_sp <- x / D
  
  # Calculate species count given transects with positive visuals
  count_sp[a, t, , 1] <- presence_all*(gamma_sp*abundance_all[a, t]*exp(nuS[a, t, y]))
  count_sp[a, t, , 2] <- presence_mature*(gamma_sp*abundance_mature[a, t]*exp(nuS[a, t, y]))
  
  return(count_sp)
}