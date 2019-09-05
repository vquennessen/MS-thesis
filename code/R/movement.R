movement <- function(t, cr, N, A, AMP) {
  
  # First area
  N[, 1, t, cr] <- (1 - AMP)*N[, 1, t, cr] + AMP*N[, 2, t, cr]
  
  # Intermediate areas
  for (a in 2:(A-1)) {
    N[, a, t, cr] <- (1 - 2*AMP)*N[, a, t, cr] + 
      AMP*(N[, a-1, t, cr] + N[, a+1, t, cr])
  }
  
  # Last area
  N[, A, t, cr] <- (1 - AMP)*N[, A, t, cr] + AMP*N[, A-1, t, cr]
  
  return(N)
  
}