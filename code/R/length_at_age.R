#' Length at age
#'
#' @param rec_age 
#' @param max_age 
#' @param L1f 
#' @param L2f 
#' @param Kf 
#' @param a1f 
#' @param a2f 
#' @param all_ages 
#'
#' @return
#' @export
#'
#' @examples
length_at_age = function(rec_age, max_age, L1f, L2f, Kf, a1f, a2f, all_ages = F) {
  
  if (all_ages == T) { 
    ages <- 1:max_age } else {
      ages <- rec_age:max_age
    }
    
  L_inf <- L1f + (L2f - L1f)/(1 - exp(-1*Kf*(a2f - a1f)))
  L <- L_inf + (L1f - L_inf)*exp(-1*Kf*(ages - a1f))
  
  return(L)
  
}




