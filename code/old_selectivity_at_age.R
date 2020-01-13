old_selectivity_at_age <- function(L, fleets, alpha, beta, start, F_fin, L50_up, 
                               L50_down, cf, switch, full) {
  
  n <- length(L)
  f <- length(fleets) # number of fleets
  selectivity <- array(rep(0, f*n), c(f, n))
  
  for (i in 1:f) {
    
    # Based on Babcock & MacCall (2011): Eq. (8)
    upcurve <- (1)/(1 + exp(-1*alpha[i]*(L - L50_up[i])))
    
    # define selectivity as asymptotic or dome-shaped for each fleet
    if (switch[i] != 0) {
      
      # Based on Babcock & MacCall (2011): Eq. (9) for dome-shaped selectivity
      downcurve <- 1 - (1 - F_fin[i])/(1 + exp(-1*beta[i]*(L - L50_down[i])))
      
    } else {
      # if asymptotic selectivity
      downcurve <- 1
    }
    
    for (j in 1:n) {
      
      # if asymptotic, or length < switch length, selectivity = upcurve
      if (switch[i] == 0 | L[j] <= switch[i]) {
        selectivity[i, j] <- upcurve[j]
        
        # otherwise, if switch length < length <= full length and dome-shaped 
        # selectivity, selectivity = 1
      } else if (L[j] > switch[i] & L[j] <= full[i] & switch[i] > 0) {
        selectivity[i, j] <- 1
        
        # finally, if length => full length and dome-shaped selectivity, 
        # selectivity = downcurve
      } else if (L[j] >= full[i] & switch[i] > 0) {
        selectivity[i, j] <- downcurve[j]
      }
      
    # multiply each row by the fraction of the fishery caught in that fleet  
    selectivity[i, ] <- cf[i]*selectivity[i, ]
      
    }
    
  }

  # Selectivity at Age
  # Based on Babcock & MacCall (2011): Eq. (7)  
  # Dimensions = 1 * age
  selectivity_at_age <- colSums(selectivity)
  
  return(selectivity_at_age)

}
