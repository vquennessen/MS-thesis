#' Title
#'
#' @param t number of timesteps, a number ( > 1)
#' @param sigma_R standard deviation of recruitment, a number ( > 0)
#' @param rho_R recruitment autocorrelation, a number (0, 1)
#'
#' @return epsilon, a numeric vector of recruitment error for each timestep
#' @export
#'
#' @examples

epsilon <- function (A, timeT, CR, NM, nuR, rho_R) {
  
  # Error term for recruitment
  # Based on Babcock & MacCall (2011): Eq. (4)
  
  # initialize epsilon vector
  # Dimensions = area * timeT * CR * M values (3)
  Eps <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM))
  
  # eps[, 1, ]
  Eps[, 1, , ] <- nuR[, 1, , ]*sqrt(1 + rho_R^2)
  
  # fill in rest of epsilon vector
  for (a in 1:A) {
    for (t in 2:timeT) {
      for (cr in 1:CR) {
        for (nm in 1:NM) {
        Eps[a, t, cr, nm] <- rho_R*Eps[a, t-1, cr, nm] + nuR[a, t, cr, nm]*sqrt(1 + rho_R^2)
        }
      }
    }
  }
  
  return(Eps)
  
}