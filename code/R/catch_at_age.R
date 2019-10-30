catch_at_age <- function(a, t, cr, nm, FM, nat_mortality, N, A, Fb, E, catch) {
  
    #' Based on the continuous (Baranov) formulation of of selectivity, i.e. that
    #' fishing mortality (FM) and natural mortality (M) are proportional to fish 
    #' abundance (N) and act simultaneously  and uniformly throughout the year, 
    #' i.e. dN/dt = -(M+F)*N
    
    coeff <- FM[ , a, t, cr, nm]/(nat_mortality[nm] + FM[ , a, t, cr, nm])
    
    catch[ , a, t, cr, nm] <- coeff * N[ , a, t, cr, nm] * exp(-1*nat_mortality[nm] 
                                                               - FM[ , a, t, cr])
  
  return(catch[, a, t, cr, nm])
  
}