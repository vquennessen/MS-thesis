#' Selectivity at age
#'
#' @param L length at age, a numeric vector
#' @param fleets fleet names, a string of character variables
#' @param alpha slope of upcurve per fleet, a numeric vector
#' @param beta slope of downcurve per fleet, a numeric vector
#' @param start length at initial vulnerability per fleet, a numeric vector
#' @param F_fin final selectivities per fleet, if dome-shaped, a numeric vector
#' @param L50_up L50 for upcurve per fleet, a numeric vector
#' @param L50_down L50 for downcurve per fleet, a numeric vector
#' @param cf fraction of fishery caught per fleet, a numeric vector
#' @param switch lengths where selectivity switches from upcurve to 1 per fleet, 
#' a numeric vector
#' @param full lengths where selectivity switches from 1 to downcurve per fleet,
#' a numeric vector
#'
#' @return
#' @export
#'
#' @examples
#' fleets <- c('sport', 'hook', 'trawl')
#' alpha <- c(0.45, 0.369, 0.426)          
#' beta <- c(1.28, 0.419, 0)              
#' start <- c(20, 20, 20)                  
#' F_fin <- c(0.265, 0.29, 0)              
#' L50_up <- c(25.2, 33.2, 46)             
#' L50_down <- c(36.7, 46.4, 50)           
#' cf <- c(0.71, 0.28, 0.01)               
#' switch <- c(29, 29, 0)                  
#' full <- c(31, 39, 0)                    
#' selectivity_at_age(L, fleets, alpha, beta, start, F_fin, L50_up, L50_down,
#' cf, switch)

selectivity_at_age <- function(L, fleets, alpha, beta, start, F_fin, L50_up, 
                               L50_down, cf, switch, full) {
  
  n <- length(L)
  f <- length(fleets) # number of fleets
  selectivity <- array(rep(0, f*n), c(f, n))
  
  for (i in 1:f) {
    
    # Based on Babcock & MacCall (2011): Eq. (8)
    upcurve <- (1)/(1 + exp(-1*alpha[i]*(L - L50_up[i])))
    
    # define selectivity as asymptotic or dome-shaped for each fleet
    if (switch[i] > 0.0001) {
      
      # Based on Babcock & MacCall (2011): Eq. (9) for dome-shaped selectivity
      downcurve <- 1 - (1 - F_fin[i])/(1 + exp(-1*beta[i]*(L - L50_down[i])))
      
    } else {
      # if asymptotic selectivity
      downcurve <- 1
    }
    
    for (j in 1:n) {
      
      # if asymptotic, or length < switch length, selectivity = upcurve
      if (switch[i] == 0 | L[j] <= switch[i]) {
        selectivity[i, j] <- upcurve[j]
        
        # otherwise, if switch length < length <= full length and dome-shaped 
        # selectivity, selectivity = 1
      } else if (L[j] > switch[i] & L[j] <= full[i] & switch[i] > 0) {
        selectivity[i, j] <- 1
        
        # finally, if length => full length and dome-shaped selectivity, 
        # selectivity = downcurve
      } else if (L[j] >= full[i] & switch[i] > 0) {
        selectivity[i, j] <- downcurve[j]
      }
      
    # multiply each row by the fraction of the fishery caught in that fleet  
    selectivity[i, ] <- cf[i]*selectivity[i, ]
      
    }
    
  }

  # Selectivity at Age
  # Based on Babcock & MacCall (2011): Eq. (7)  
  # Dimensions = 1 * age
  selectivity_at_age <- colSums(selectivity)
  
  return(selectivity_at_age)

}
