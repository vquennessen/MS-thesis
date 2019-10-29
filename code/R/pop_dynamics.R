pop_dynamics <- function(a, t, cr, rec_age, max_age, n, SSB, N, W, Mat, A, R0, 
                         h, B0, Eps, sigma_R, Fb, E, S, M, FM, m, abundance_all, 
                         abundance_mature, biomass, fishing) {
  
  # Calculate fishing mortality
  if (fishing == T) {
    FM[, a, t, cr] <- fishing_mortality(a, t, cr, FM, A, Fb, E, S)
  } else { FM[, a, t, cr] <- 0 }
  
  ##### Step population foward in time
  
  # Calculate recruitment and add recruits to population
  N[1, a, t, cr] <- recruitment(a, t, cr, SSB, A, R0, h, B0, Eps, sigma_R)
  
  # Ages rec_age + 1 to max_age - 1
  for (i in 2:(n - 1)) {
    N[i, a, t, cr] <- N[i - 1, a, t - 1, cr] * exp(-1 * (FM[i - 1, a, t - 1, cr] + M))
  }
  
  # Final age bin
  N[n, a, t, cr] <- N[n - 1, a, t - 1, cr] * exp(-1 * (FM[n - 1, a, t - 1, cr] + M)) + 
    N[n, a, t - 1, cr] * exp(-1 * (FM[n, a, t - 1, cr] + M)) 
  
  # Calculate abundance of all fish
  abundance_all[a, t, cr] <- sum(N[, a, t, cr])
  
  # Calculate abundance of mature fish
  abundance_mature[a, t, cr] <- sum(N[m:n, a, t, cr])
  
  # Calculate biomass of all fish
  biomass[a, t, cr] <- sum(N[, a, t, cr] * W)
  
  # Calculate spawning stock biomass
  SSB[a, t, cr] <- sum(N[, a, t, cr]*W*Mat)
  
  output <- list(SSB, FM, N, abundance_all, abundance_mature, biomass)
  
  return(output)
  
}