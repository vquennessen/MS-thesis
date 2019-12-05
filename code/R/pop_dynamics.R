pop_dynamics <- function(a, t, cr, nm, rec_age, max_age, n, SSB, N, W, Mat, A, 
                         R0, h, B0, Eps, sigma_R, Fb, E, S, NM, FM, m, 
                         abundance_all, abundance_mature, biomass, fishing, 
                         nat_mortality, recruitment_mode) {
  
  # Calculate fishing mortality
  if (fishing == T) {
    FM[, a, t, cr, nm] <- fishing_mortality(a, t, cr, nm, FM, A, Fb, E, S)
  } else { FM[, a, t, cr, nm] <- 0 }
  
  ##### Step population foward in time
  
  # Calculate recruitment and add recruits to population
  N[1, a, t, cr, nm] <- recruitment(a, t, cr, nm, SSB, A, R0, h, B0, Eps, 
                                    sigma_R, rec_age, recruitment_mode)
  
  # Ages rec_age + 1 to max_age - 1
  for (i in 2:(n - 1)) {
    N[i, a, t, cr, nm] <- N[i - 1, a, t - 1, cr, nm] * exp(-1 * (FM[i - 1, a, t - 1, cr, nm] + nat_mortality[nm]))
  }
  
  # Final age bin
  N[n, a, t, cr, nm] <- N[n - 1, a, t - 1, cr, nm] * exp(-1 * (FM[n - 1, a, t - 1, cr, nm] + nat_mortality[nm])) + 
    N[n, a, t - 1, cr, nm] * exp(-1 * (FM[n, a, t - 1, cr, nm] + nat_mortality[nm])) 
  
  # Calculate abundance of all fish
  abundance_all[a, t, cr, nm] <- sum(N[, a, t, cr, nm])
  
  # Calculate abundance of mature fish
  abundance_mature[a, t, cr, nm] <- sum(N[m:max_age - 1, a, t, cr, nm])
  
  # Calculate biomass of all fish
  biomass[a, t, cr, nm] <- sum(N[, a, t, cr, nm] * W)
  
  # Calculate spawning stock biomass
  SSB[a, t, cr, nm] <- sum(N[, a, t, cr, nm]*W*Mat)
  
  output <- list(FM[, a, t, cr, nm], N[, a, t, cr, nm], 
                 abundance_all[a, t, cr, nm], abundance_mature[a, t, cr, nm], 
                 biomass[a, t, cr, nm], SSB[a, t, cr, nm])
  
  return(output)
  
}