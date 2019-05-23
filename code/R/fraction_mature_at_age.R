#' Fraction mature at age
#' 
#' @param max_age  maximum age, a number (> 0)
#' @param k_mat    slope of maturity curve, a number
#' @param L        length at age, a vector of numbers
#' @param L50      length at 50% maturity, a number
#' 
#' @return M       fraction mature at age, a vector of numbers
#' 
#' @example 
#' L <- length_at_age(35, 32.21, 47.95, 0.2022, 5, 15)
#' fraction_mat_at_age(35, -0.4103, L, 39.53)

fraction_mature_at_age = function(n, k_mat, L, L50) {
  
  M <- array(rep(0, n), c(1, n))
  M <- (1)/(1 + exp(k_mat*(L - L50)))
  
  return(M)
  
}