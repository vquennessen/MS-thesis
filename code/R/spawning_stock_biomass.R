#' Spawning Stock Biomass
#'
#' @param N population age structure, a numeric vector
#' @param W weight at age, a numeric vector
#' @param M fraction mature at age, a numeric vector
#'
#' @return SSB spawning stock biomass for the population
#' @export
#'
#' @examples
#' par <- parameters("black rockfish")
#' max_age <- par[[1]]
#' a1f <- par[[8]]; L1f <- par[[9]]; a2f <- par[[10]]; L2f <- par[[11]]; 
#' Kf <- par[[12]]; af <- par[[4]]; bf <- par[[5]]
#' L50 <- par[[18]]; k_mat <- par[[19]]  
#' 
#' N <- array(rep(100, max_age), c(1, max_age))
#' L <- length_at_age(max_age, L1f, L2f, Kf, a1f, a2f)
#' W <- weight_at_age(L, af, bf)
#' Mat <- fraction_mature_at_age(max_age, k_mat, L, L50)
#' 
#' spawning_stock_biomass(N, W, M)at

spawning_stock_biomass <- function(a, t, cr, N, W, Mat) {
  
  # Dimensions = area * time * CR
  SSB[a, t, cr] <- sum(N[, a, t - rec_age, cr]*W*Mat)
  
  return(SSB[a, t, cr])
  
}