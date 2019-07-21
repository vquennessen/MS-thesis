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
#' par <- parameters("black rockfish")
#' max_age <- par[[1]]
#' a1f <- par[[8]]; L1f <- par[[9]]; a2f <- par[[10]]; L2f <- par[[11]]; 
#' Kf <- par[[12]]; af <- par[[4]]; bf <- par[[5]]
#' L50 <- par[[18]]; k_mat <- par[[19]]
#' R0 <- par[[21]]; h <- par[[22]]; sigma_R <- par[[24]]; B0 <- par[[32]] 
#' 
#' A <- 5
#' 
#' epsilon <- epsilon(1, sigma_R, rho_R)
#' N <- array(rep(100, max_age), c(1, max_age))
#' L <- length_at_age(max_age, L1f, L2f, Kf, a1f, a2f)
#' W <- weight_at_age(L, af, bf)
#' M <- fraction_mature_at_age(max_age, k_mat, L, L50)
#' B <- spawning_stock_biomass(N, W, M)

recruitment = function(a, t, cr, SSB, A, R0, h, B0, Eps, sigma_R) {
  
  # Recruitment
  # Based on Babcock & MacCall (2011): Eq. (3)
  # Dimensions = 1 * 1
  
  adjR0 <- R0 / A
  
  R1 <- (0.8 * adjR0 * h * SSB[a, t-1, cr]) / (0.2 * B0 * (1 - h) + (h - 0.2) * SSB[a, t-1, cr]) 
  R <- R1 * (exp(Eps - sigma_R^2 / 2))

  return(R)
  
}