# Selectivity at age per fleet

selectivity <- function(L, fleets, alpha, beta, start, F_fin, L50_up, L50_down) {
  
  f <- length(fleets) # number of fleets
  
  # Based on Babcock & MacCall (2011): Eq. (8)
  upcurve <- (1)/(1 + exp(-1*alpha[1]*(L - L50_up)))

  
  if (switch > 0.0001) {
    
    # Based on Babcock & MacCall (2011): Eq. (9) for dome-shaped selectivity
    downcurve <- 1 - (1 - F_fin)/(1 + exp(-1*beta[1]*(L - L50_down)))
  
  } else {
    # if asymptotic selectivity
    downcurve <- 1
  }

  # TODO: What is switch?!?
  
  for (i in 1:f) {
    
  }
  
  
  # Selectivity at Age
  # Based on Babcock & MacCall (2011): Eq. (7)
  Sp <- array(rep(0, n*f), c(n, f))

}