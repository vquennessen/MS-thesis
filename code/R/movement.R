movement <- function(t, cr, nm, N, A, AMP) {
  
  # First area to second area
  N[, 1, t, cr, nm] <- (1 - AMP)*N[, 1, t, cr, nm] + AMP*N[, 2, t, cr, nm]
  
  # Intermediate areas to adjacent areas
  for (a in 2:(A-1)) {
    N[, a, t, cr, nm] <- (1 - 2*AMP)*N[, a, t, cr, nm] + 
      AMP*(N[, a - 1, t, cr, nm] + N[, a + 1, t, cr, nm])
  }
  
  # Last area to next to last area
  N[, A, t, cr, nm] <- (1 - AMP)*N[, A, t, cr, nm] + AMP*N[, A - 1, t, cr, nm]
  
  return(N)
  
}