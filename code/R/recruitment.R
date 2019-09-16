#' Recruitment
#'
#' @param SSB spawning stock biomass, a number ( > 0)
#' @param A number of areas, a number ( > 0)
#' @param R0 unfished recruitment, a number ( > 0)
#' @param h steepness, a number ( > 0)
#' @param B0 unfished spawning stock biomass, a number ( > 0)
#' @param epsilon an error term, a number 
#' @param sigma_R the standard deviation of a normally distributed random 
#' variable, a number
#'
#' @return R recruitment
#' @export
#'
#' @examples

recruitment = function(a, t, cr, SSB, A, R0, h, B0, Eps, sigma_R) {
  
  # Recruitment
  # Based on Babcock & MacCall (2011): Eq. (3)
  # Dimensions = 1 * 1
  
  adjR0 <- R0 / A
  
  R1 <- (0.8 * adjR0 * h * SSB[a, t-1, cr]) / (0.2 * B0 * (1 - h) + (h - 0.2) * SSB[a, t-1, cr]) 
  R <- R1 * (exp(Eps - sigma_R^2 / 2))

  return(R)
  
}