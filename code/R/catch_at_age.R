catch_at_age <- function(a, t, cr, N, FM, catch_formulation, fishing) {
  
  if (catch == 'continuous') {
    
    #' Based on the continuous (Baranov) formulation of of selectivity, i.e. that
    #' fishing mortality (FM) and natural mortality (M) are proportional to fish 
    #' abundance (N) and act simultaneously  and uniformly throughout the year, 
    #' i.e. dN/dt = -(M+F)*N
    
    cf <- FM[ , a, t, cr]/(M + FM[ , a, t, cr])
    
    catch[ , a, t, cr] <- cf * N[ , a, t, cr]*exp(-1*M - FM[ , a, t, cr])
    
  } else if (catch == 'discrete') {
    
    #' Based on the discrete formulation of selectivity, which assumes that 
    #' fishing occurs in a very short pulse and no natural mortality occurs
    #' during the fishing season, which could be at the start, end, or any 
    #' other part of the year, determined by the variable fishing, which is
    #' 0 for at the beginning of the year, 1 for the end, and on the interval
    #' (0, 1) for anywhere in between
    
    u <- 1 - exp(FM[, a, t, cr]) 
    
    catch[, a, t, cr] <- N[, a, t, cr]*S*u*exp(-1*(fishing - t)*M)
    
  }
  
  return(catch)
  
}