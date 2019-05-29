sampling <- function(a, t, r, D, abundance, transects) {
  
  # Calculate delta
  # Based on Babcock & MacCall (2011): Eq. (13)
  delta <- r / D
  
  # Calculate probability of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  p <-  delta * abundance[a, t]
  
  # Determine if species is seen at least once
  # Dimensions = 1 * transects
  presence <- rbinom(transects, 1, p)
  
  # Calculate species count given transects with positive visuals
  
  
}