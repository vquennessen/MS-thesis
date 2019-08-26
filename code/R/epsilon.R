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
#' t <- 50 
#' par <- parameters("black rockfish")
#' sigma_R <- par[[24]]
#' rho_R <- par[[25]]
#' epsilon(t, sigma_R, rho_R)

epsilon <- function (A, time, CR, nuR, sigma_R, rho_R) {
  
  # Error term for recruitment
  # Based on Babcock & MacCall (2011): Eq. (4)
  
  # initialize epsilon vector
  # Dimensions = area * time * CR
  Eps <- array(rep(0, A*time*CR), c(A, time, CR))
  
  # eps[, 1, ]
  Eps[, 1, ] <- nuR[, 1, ]*sqrt(1 + rho_R^2)
  
  # fill in rest of epsilon vector
  for (a in 1:A) {
    for (t in 2:time) {
      for (cr in 1:CR) {
        Eps[a, t, cr] <- rho_R*Eps[a, t-1, cr] + nuR[a, t, cr]*sqrt(1 + rho_R^2)
      }
    }
  }
  
  return(Eps)
  
}