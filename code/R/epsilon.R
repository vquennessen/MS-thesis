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

epsilon <- function (A, t, CR, sigma_R, rho_R) {
  
  # Error term for recruitment
  # Based on Babcock & MacCall (2011): Eq. (4)
  
  # initialize epsilon vector
  # Dimensions = area * time * CR
  epsilon <- array(rep(0, A*t*CR), c(A, t, CR))
  
  # sample from normal distribution with mean 0 and std sigma_R
  nu <- array(rnorm(A*t*CR, 0, sigma_R), c(A, t, CR))
  
  # epsilon[, 1, ]
  epsilon[, 1, ] <- nu[, 1, ]*sqrt(1 + rho_R^2)
  
  # fill in rest of epsilon vector
  for (a in 1:A) {
    for (t in 2:t) {
      for (y in 1:CR) {
        epsilon[a, t, y ] <- rho_R*epsilon[a, t-1, y] + nu[a, t, y]*sqrt(1 + rho_R^2)
      }
    }
  }
  
  return(epsilon)
  
}