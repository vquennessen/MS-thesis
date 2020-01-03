#' Length at age
#' Based on Babcock & MacCall (2011): Eq. (10)
#'
#' @param max_age   maximum age, a number (> 0)
#' @param L1f       length at age 1, female, a number (> 0)
#' @param L2f       length at age 2, female, a number (> L1f)
#' @param Kf        von Bertalanffy growth rate, female, a number (0, 1)
#' @param a1f       age 1, female, a number (> 0)
#' @param a2f       age 2, female, a number (> a1f)

length_at_age = function(rec_age, max_age, L1f, L2f, Kf, a1f, a2f, all_ages) {
  
  ages <- ifelse(all_ages == T, 1:max_age, rec_age:max_age)
  
  L_inf <- L1f + (L2f - L1f)/(1 - exp(-1*Kf*(a2f - a1f)))
  L <- L_inf + (L1f - L_inf)*exp(-1*Kf*(ages - a1f))
  
  return(L)
  
}

