catch_at_age <- function(a, t, cr, nm, FM, Nat_mortality, N, A, Fb, E, Catch) {
  
    #' Based on the continuous (Baranov) formulation of of selectivity, i.e. that
    #' fishing mortality (FM) and natural mortality (M) are proportional to fish 
    #' abundance (N) and act simultaneously  and uniformly throughout the year, 
    #' i.e. dN/dt = -(M+F)*N
    
    coeff <- FM[ , a, t, cr, nm]/(Nat_mortality[nm] + FM[ , a, t, cr, nm])
    
    Catch[ , a, t, cr, nm] <- coeff * N[ , a, t, cr, nm] * exp(-1*Nat_mortality[nm] 
                                                               - FM[ , a, t, cr, nm])
  
  return(Catch[, a, t, cr, nm])
  
}