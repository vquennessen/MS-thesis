sampling <- function(a, t, r, D, abundance, transects, x, count_sp) {
  
  # Calculate delta
  # Based on Babcock & MacCall (2011): Eq. (13)
  delta <- r / D
  
  # Calculate probability of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  p <-  delta * abundance[a, t]
  
  # Determine if species is seen at least once
  # Dimensions = 1 * transects
  presence <- rbinom(transects, 1, p)
  
  # Based on Babcock & MacCall (2011): Eq. (16)
  gamma_sp <- x / D
  
  # Calculate standard deviation of normal variable
  # Based on Babcock & MacCall (2011): Eq. (15)
  sigma_sp <- sqrt(log(1 + (sp/x)^2))
  
  # TODO generate all rnorm values, then call on them for next calculation
  
  # Calculate species count given transects with positive visuals
  count_sp[a, t] <- presence*(gamma_sp*abundance[a, t]*exp(rnorm(1, 0, sigma_sp)))
  
  return(count_sp)
}