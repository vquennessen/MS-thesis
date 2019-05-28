sampling <- function(a, t, r, D, p, abundance) {
  
  # calculate delta
  # Based on Babcock & MacCall (2011): Eq. (13)
  delta <- r / D
  
  # calculate probability of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  p[a, t] <-  delta * abundance[a, t]
  
  # determine if species is seen at least once
  observations <- 
  trials <- 10
  
  presence[b] <- 
  
  # calculate species count given transects with positive visuals
  
  
}