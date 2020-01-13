epsilon <- function (A, timeT, CR, NM, NuR, Rho_R) {
  
  # Error term for recruitment
  # Based on Babcock & MacCall (2011): Eq. (4)
  
  # initialize epsilon vector
  # Dimensions = area * timeT * CR * M values (3)
  Eps <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM))
  
  # eps[, 1, ]
  Eps[, 1, , ] <- NuR[, 1, , ]*sqrt(1 + Rho_R^2)
  
  # fill in rest of epsilon vector
  for (a in 1:A) {
    for (t in 2:TimeT) {
      for (cr in 1:CR) {
        for (nm in 1:NM) {
        Eps[a, t, cr, nm] <- Rho_R*Eps[a, t-1, cr, nm] + NuR[a, t, cr, nm]*sqrt(1 + Rho_R^2)
        }
      }
    }
  }
  
  return(Eps)
  
}