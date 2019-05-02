#' Weight at Age
#'
#' @param L       length at age, a vector of numbers (> 0)
#' @param af      growth parameter, a number (0 - 1)
#' @param bf      growth parameter, a number (~ 3)
#'
#' @return Weight at age vector
#'
#' @examples 
#' L <- length_at_age(35, 32.21, 47.95, 0.2022, 5, 15)
#' weight_at_age(L, 1.68e-5, 3)

weight_at_age = function(L, af, bf) {
  
  # Weight at age
  # Based on Babcock & MacCall (2011): Eq. (11)
  W <- af*L^bf
  
  return(W)
  
}