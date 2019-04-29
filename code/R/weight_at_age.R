#' Weight at Age
#'
#' @param max_age maximum age, a number (> 1)
#' @param a1f     age 1, a number (> 0)
#' @param a2f     age 2, a number (> age 1)
#' @param L1f     length at age 1, a number (> 0)
#' @param L2f     length at age 2, a number (> length at age 1)
#' @param Kf      von Bertalanffy growth rate, a number (0 - 1)
#' @param af      growth parameter, a number (0 - 1)
#' @param bf      growth parameter, a number (~ 3)
#'
#' @return Weight at age vector
#' @export
#'
#' @examples
#' weight_at_age(35, 5, 15, 32.21, 47.95, 0.2022, 1.68e-5, 3)

weight_at_age = function(max_age, a1f, a2f, L1f, L2f, Kf, af, bf) {
  
  # Length at age
  # Based on Babcock & MacCall (2011): Eq. (10)
  age <- 1:max_age
  L_inf <- L1f + (L2f - L1f)/(1 - exp(-1*Kf*(a2f - a1f)))
  L <- L_inf + (L1f - L_inf)*exp(-1*Kf*(age - a1f))
  
  # Weight at age
  # Based on Babcock & MacCall (2011): Eq. (11)
  W <- af*L^bf
  
  return(W)
  
}