#' Title
#'
#' @param A number of areas, a number
#' @param Fb fishing mortality rate before reserve implementation, a number
#' @param E fishing effort in each area before reserve implementation, a number
#' @param S selectivity at age, a numeric vector
#'
#' @return fishing_M, fishing mortality over age and area, a 2D numeric array
#' @export
#'
#' @examples
#' 
#' par <- parameters("black rockfish")
#' max_age <- par[[1]]; Fb <- par[[28]]
#' a1f <- par[[8]]; L1f <- par[[9]]; a2f <- par[[10]]; L2f <- par[[11]]; 
#' Kf <- par[[12]]; a1m <- par[[13]]; L1m <- par[[14]];
#' a2m <- par[[15]]; L2m <- par[[16]]; Km <- par[[17]]  
#' L50 <- par[[18]]; fleets <- par[[32]]; alpha <- par[[33]]
#' beta <- par[[34]]; start <- par[[35]]; F_fin <- par[[36]]
#' L50_up <- par[[37]]; L50_down <- par[[38]]; cf <- par[[39]]
#' switch <- par[[40]]; full <- par[[41]]; A <- 5; E <- 0.10
#' 
#' L <- length_at_age(max_age, L1f, L2f, Kf, a1f, a2f)
#' 
#' S <- selectivity_at_age(L, fleets, alpha, beta, start, F_fin, L50_up, 
#' L50_down, cf, switch, full)
#' 
#' fishing_mortality(A, Fb, E, S)

fishing_mortality <- function(a, t, cr, FM, A, Fb, E, S) {

  # Catchability
  # Based on Babcock & MacCall (2011): Eq. (6)
  catchability <- (A*Fb)/(sum(E[, t, cr]))
  
  # Selectivity as a matrix
  # dimensions = age * 1
  selectivity <- array(S, c(length(S), 1)) 
  
  # Effort as a matrix
  # Dimensions = 1 * area
  effort <- E[a, t, cr]

  # Fishing mortality
  # Based on Babcock & MacCall (2011): Eq. (5)
  # Dimensions = age * area * time * CR
  FM[, a, t, cr] <- catchability * selectivity * effort
  
  return(FM[, a, t, cr])
  
}