#' Based on the continuous (Baranov) formulation of of selectivity, i.e. that
#' fishing mortality (FM) and natural mortality (M) are proportional to fish 
#' abundance (N) and act simultaneously  and uniformly throughout the year, 
#' i.e. dN/dt = -(M+F)*N

catch_at_age <- function(a, t, cr, N, FM) {
  
  cf <- FM[ , a, t, cr]/(M + FM[ , a, t, cr])
  
  catch[ , a, t, cr] <- cf * N[ , a, t, cr]*exp(-1*M - FM[ , a, t, cr])
  
  return(catch)
  
}