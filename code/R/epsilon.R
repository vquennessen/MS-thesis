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

epsilon <- function (t, sigma_R, rho_R) {
  
  # Error term for recruitment
  # Based on Babcock & MacCall (2011): Eq. (4)
  
  # initialize epsilon vector
  epsilon <- array(rep(0, t), t)
  
  # sample from normal distribution with mean 0 and std sigma_R
  nu <- rnorm(t, 0, sigma_R)
  
  # epsilon[1]
  epsilon[1] <- nu[1]*sqrt(1 + rho_R^2)
  
  # fill in rest of epsilon vector
  for (i in 2:t) {
    epsilon[i] <- rho_R*epsilon[i-1] + nu[i]*sqrt(1 + rho_R^2)
  }
  
  return(epsilon)
  
}