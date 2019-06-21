#' Population Dynamics, step forward one year in time
#'
#' @param a area, a number (1, A)
#' @param t time step, a number (1, time)
#' @param SSB spawning stock biomass, a 2D numeric array (dimensions area x time)
#' @param N age-structured population, a 3D numeric array 
#' (dimensions age x area x time)
#' @param W weight at age, a numeric vector (length age)
#' @param Mat fraction mature at age, a numeric vector (length age)
#' @param A total number of areas, a number ( > 0)
#' @param R0 unfished recruitment, a number ( > 0)
#' @param h steepness parameter, a number
#' @param B0 unfished spawning stock biomass, a number ( > 0)
#' @param e epsilon, a numeric vector (length time)
#' @param sigma_R recruitment standard deviation, a number ( > 0)
#' @param Fb fully selected fishing mortality rate on the whole population 
#' before reserve implementation, a number (0, 1)
#' @param E nominal fishing effort in each area before reserve implementation,
#' a number (0, 1)
#' @param S selectivity at age, a numeric vector (length age)
#' @param M natural mortality, a number (0, 1)
#'
#' @return B, R, FM, N - updated spawning stock biomass, recruitment, fishing
#' mortality, and age-structured population size, after 1 time step
#' @export
#'
#' @examples

pop_dynamics <- function(a, t, y, rec_age, max_age, n, SSB, N, W, Mat, A, R0, h, 
                         B0, e, sigma_R, Fb, E, S, M) {
  
  # Calculate spawning stock biomass
  SSB[a, t] <- spawning_stock_biomass(N[, a, t-rec_age], W, Mat)
  
  # Calculate recruitment and add recruits to population
  N[1, a, t] <- recruitment(SSB[a, t-1], A, R0, h, B0, e[a, t, y], sigma_R)
  
  # Calculate fishing mortality
  FM[, , t] <- fishing_mortality(A, Fb, E, S)
  
  # Step population foward in time
  for (i in 2:(n - 1)) {
    N[i, a, t] <- N[i - 1, a, t - 1] * exp(-1 * (FM[i - 1, a, t - 1] + M))
  }
  
  N[max_age - 1, a, t] <- N[(max_age - 2), a, t - 1] * 
    exp(-1 * (FM[(max_age - 2), a, t - 1] + M)) + 
    N[max_age - 1, a, t - 1] * exp(-1 * (FM[max_age - 1, a, t - 1] + M)) 
  
  output <- list(SSB, R, FM, N)
  
  return(output)
  
}