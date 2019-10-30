fishing_mortality <- function(a, t, cr, nm, FM, A, Fb, E, S) {
  
  # Catchability (Vulnerability to fishing gear)
  # Based on Babcock & MacCall (2011): Eq. (6)
  vulnerability <- (A*Fb)/(sum(E[, 1, 1, 1]))
  
  # Selectivity as a matrix
  # dimensions = age * 1
  selectivity <- array(S, c(length(S), 1)) 
  
  # Effort as a matrix
  # Dimensions = area * time * CR
  effort <- E[a, t, cr, nm]
  
  # Fishing mortality
  # Based on Babcock & MacCall (2011): Eq. (5)
  # Dimensions = age * area * time * CR
  FM[, a, t, cr, nm] <- vulnerability * selectivity * effort
  
  return(FM[, a, t, cr, nm])
  
}