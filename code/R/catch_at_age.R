#' Catch at age
#'
#' @param a area, a number (1, A)
#' @param t time, a number (1, timeT)
#' @param cr control rule, a number (1, CR)
#' @param FM fishing mortality, a 1D numeric array (length = age)
#' @param M natural mortality, a number (0, 1)
#' @param N population size, a 4D numeric array 
#' (dimensions = age*area*time*control rule)
#' @param A total number of areas, a number (> 1)
#' @param Fb initial fishing mortality that leads to depletion, a number (0, 1)
#' @param E fishing effort, a 3D numeric array 
#' (dimensions = area*time*control rule)
#' @param catch catch at age, a 4D numeric array 
#' (dimensions = age*area*time*control rule)
#' @param catch_form catch formulation, character, either 'continuous' or 
#' 'discrete'
#' @param season the season over which fishing occurs, 
#' only applies if catch_form = 'discrete', a number (0, 1)
#'
#' @return returns the updated catch at age for all ages in area *a*, at time 
#' *t*, under control rule *cr*
#' @export
#'
#' @examples
catch_at_age <- function(a, t, cr, FM, M, N, A, Fb, E, catch, catch_form, 
                         season) {
  
  if (catch_form == 'continuous') {  
    
    #' Based on the continuous (Baranov) formulation of of selectivity, i.e. that
    #' fishing mortality (FM) and natural mortality (M) are proportional to fish 
    #' abundance (N) and act simultaneously  and uniformly throughout the year, 
    #' i.e. dN/dt = -(M+F)*N
    
    coeff <- FM[ , a, t, cr]/(M + FM[ , a, t, cr])
    
    catch[ , a, t, cr] <- coeff * N[ , a, t, cr] * exp(-1*(M + FM[ , a, t, cr]))
    
  } else if (catch_form == 'discrete') {
    
    #' Based on the discrete formulation of selectivity, which assumes that
    #' fishing occurs in a very short pulse and no natural mortality occurs
    #' during the fishing season, which could be at the start, end, or any
    #' other part of the year, determined by the variable fishing, which is
    #' 0 for at the beginning of the year, 1 for the end, and on the interval
    #' (0, 1) for anywhere in between
    
    vulnerability <- vulnerability_to_gear(a, t, cr, A, Fb, E)
    
    u <- 1 - exp(-1*vulnerability*E[a, t, cr])
    
    catch[, a, t, cr] <- N[, a, t, cr]*S*u*exp(-1*(season - t)*M)
    
  }
  
  return(catch[, a, t, cr])
  
}