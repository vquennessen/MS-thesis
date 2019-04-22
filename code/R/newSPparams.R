#' Loads life history characteristic parameters for new species
#' 
#' @param max_age maximum age, a number (> 0)
#' @param M natural mortality, a number (0, 1)
#' @param rec_age age at recruitment, a number (> 0, < max_age)
#' @param af female weight at length parameter, a number (> 0)
#' @param bf female weight at length parameter, a number (~ 3)
#' @param am male weight at length parameter, a number (> 0)
#' @param bm male weight at length parameter, a number (~ 3)
#' @param a1f female age 1, a number (> 0)
#' @param L1f female length at age 1, a number (> 0)
#' @param a2f female age 2, a number (> a1f)
#' @param L2f female length at age 2, a number (> L1f)
#' @param Kf female growth rate, a number (0, 1)
#' @param a1m male age 1, a number (> 0)
#' @param L1m male length at age 1, a number (> 0)
#' @param a2m male age 2, a number (> a1m)
#' @param L2m male length at age 2, a number (> L1m)
#' @param Km male growth rate, a number (0, 1)
#' @param L50 length at 50% maturity, a number (> 0)
#' @param k_mat slope of maturity curve, a number
#' @param ldp larval drift proportion, a number (0, 1)
#' @param R0 unfished recruitment, a number
#' @param h steepness, a number (0, 1)
#' @param phi unfished recruits per spawner, a number
#' @param sigma_R recruitment standard deviation, a number
#' @param rho_R recruitment autocorrelation, a number
#' @param p adult movement proportion, a number (0, 1)
#' @param D depletion, a number (0, 1)
#' @param Fb fishing mortality to cause D, a number (0, 1)
#' @param r proportion of positive transects in monitoring data, a number (0, 1)
#' @param x mean of positive transects, a number
#' @param sp standard deviation of positive transects, a number
#' 
#' @return the parameters needed to run the popDyn function

newSPparams = function(max_age, M, rec_age, af, bf, am, bm, a1f, L1f, a2f, L2f, 
                       Kf, a1m, L1m, a2m, L2m, Km, L50, k_mat, ldp, R0, h, phi,
                       sigma_R, rho_R, p, D, Fb, r, x, sp) {

  }