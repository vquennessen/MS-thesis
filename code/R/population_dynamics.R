population_dynamics <- function(a, t, B, N, W, M, A, R0, h, B0, e, sigma_R, Fb, E, S) {
  
  # Calculate spawning stock biomass
  B[a, t] <- spawning_stock_biomass(N[, a, t-1] * W * M)
  
  # Calculate recruitment
  R <- recruitment(B[a, t], A, R0, h, B0, e[t], sigma_R)
  
  # Add recruits to population
  N[1, a, t] <- R
  
  # Calculate fishing mortality
  FM[, , t] <- fishing_mortality(A, Fb, E, S)
  
  # Step population foward in time
  for (i in age) {
    N[i, a, t] <- N[i-1, s, t-1] * exp(-1 * FM[i-1, a, t-1] + M)
  }
  
  return(B, R, FM, N)
  
}