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

fishing_mortality <- function(A, Fb, E, S) {
  
  # dimensions = age * area
  selectivity <- array(rep(S, A), c(length(S), A))
  
  # dimensions = area * area
  effort <- array(rep(E, A), c(A, A))
  
  # Catchability
  # Based on Babcock & MacCall (2011): Eq. (6)
  catchability <- (A*Fb)/(sum(effort[1, ]))
  
  # Fishing mortality
  # Based on Babcock & MacCall (2011): Eq. (5)
  # Dimensions = age * area
  fishing_M <- catchability * selectivity %*% effort
  
  return(fishing_M)
  
}