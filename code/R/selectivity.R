# Selectivity at age per fleet

selectivity <- function(L, fleets, alpha, beta, start, F_fin, L50_up, L50_down, 
                        switch) {
  
  n <- length(L)
  f <- length(fleets) # number of fleets
  selectivity <- array(rep(0, n*f), c(f, n))
  
  for (i in 1:f) {
    
    # Based on Babcock & MacCall (2011): Eq. (8)
    upcurve <- (1)/(1 + exp(-1*alpha[i]*(L - L50_up[i])))
    
    # define selectivity as asymptotic or dome-shaped for each fleet
    if (switch > 0.0001) {
      
      # Based on Babcock & MacCall (2011): Eq. (9) for dome-shaped selectivity
      downcurve <- 1 - (1 - F_fin[i])/(1 + exp(-1*beta[i]*(L - L50_down[i])))
      
    } else {
      # if asymptotic selectivity
      downcurve <- 1
    }
    
    for (j in 1:n) {
      if (switch == 0 | L[j] <= switch[i]) {
        selectivity[i, j] <- upcurve[j]
      } else if (L[j] > switch[i] & L[j] <= full[i] & switch[i] > 0) {
        selectivity[i, j] <- 1
      } else if (L[j] >= full[i] & switch[i] > 0) {
        selectivity[i, j] <- downcurve[j]
      }
      
    }
    
  }

  # Selectivity at Age
  # Based on Babcock & MacCall (2011): Eq. (7)

}
