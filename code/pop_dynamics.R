pop_dynamics <- function(a, t, cr, nm, Rec_age, Max_age, Num, SSB, N, W, Mat, A, 
                         R0, H, B0, Eps, Sigma_R, Fb, E, S, NM, FM, A50_mat, 
                         Abundance_all, Abundance_mature, Biomass, Fishing, 
                         Nat_mortality, Recruitment_mode) {
  
  # Calculate fishing mortality
  if (Fishing == T) {
    FM[, a, t, cr, nm] <- fishing_mortality(a, t, cr, nm, FM, A, Fb, E, S)
  } else { FM[, a, t, cr, nm] <- 0 }
  
  ##### Step population foward in time
  
  # Calculate recruitment and add recruits to population
  N[1, a, t, cr, nm] <- recruitment(a, t, cr, nm, SSB, A, R0, H, B0, Eps, 
                                    Sigma_R, Rec_age, Recruitment_mode)
  
  # Ages rec_age + 1 to max_age - 1
  for (i in 2:(Num - 1)) {
    N[i, a, t, cr, nm] <- N[i - 1, a, t - 1, cr, nm] * exp(-1 * (FM[i - 1, a, t - 1, cr, nm] + Nat_mortality[nm]))
  }
  
  # Final age bin
  N[Num, a, t, cr, nm] <- N[Num - 1, a, t - 1, cr, nm] * exp(-1 * (FM[Num - 1, a, t - 1, cr, nm] + Nat_mortality[nm])) + 
    N[Num, a, t - 1, cr, nm] * exp(-1 * (FM[Num, a, t - 1, cr, nm] + Nat_mortality[nm])) 
  
  # Calculate abundance of all fish
  Abundance_all[a, t, cr, nm] <- sum(N[, a, t, cr, nm])
  
  # Calculate abundance of mature fish
  Abundance_mature[a, t, cr, nm] <- sum(N[A50_mat:Max_age - 1, a, t, cr, nm])
  
  # Calculate biomass of all fish
  Biomass[a, t, cr, nm] <- sum(N[, a, t, cr, nm] * W)
  
  # Calculate spawning stock biomass
  SSB[a, t, cr, nm] <- sum(N[, a, t, cr, nm]*W*Mat)
  
  output <- list(FM[, a, t, cr, nm], N[, a, t, cr, nm], 
                 Abundance_all[a, t, cr, nm], Abundance_mature[a, t, cr, nm], 
                 Biomass[a, t, cr, nm], SSB[a, t, cr, nm])
  
  return(output)
  
}