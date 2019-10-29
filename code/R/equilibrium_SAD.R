equilibrium_SAD <- function(a, cr, A, rec_age, max_age, n, W, R0,
                            Mat, h, B0, Eps, sigma_R, Fb, S, M, eq_time, m, 
                            stochasticity, rho_R) {
  
  # Initialize population size array
  # Dimensions = age * 1 * time * 1
  N2 <- FM2 <- array(rep(0, n*eq_time), c(n, 1, eq_time, 1))
  
  # Initialize biomass and SSB arrays
  # Dimensions = 1 * time * 1
  biomass2 <- SSB2 <- array(rep(0, eq_time), c(1, eq_time, 1))
  abundance_all2 <- abundance_mature2 <- array(rep(0, eq_time), c(1, eq_time, 1))
  
  # Recruitment normal variable
  # Dimensions = area * timeT * CR
  if (stochasticity == T) {
    nuR2 <- array(rnorm(eq_time, 0, sigma_R), c(1, eq_time, 1))
  } else if (stochasticity == F) {
    nuR2 <- array(rep(0, 1*eq_time*1), c(1, eq_time, 1))
  }
  
  # Recruitment error
  # Dimensions = area * timeT * CR
  Eps2 <- epsilon(A = a, eq_time, CR = cr, nuR2, rho_R)
  
  # Start each age class with 10 individuals
  # Enter FM, N, abundance, and biomasses for time = 1 to rec_age
  # Dimensions = age * area * time * CR
  for (t in 1:(rec_age + 1)) {
    N2[, a, t, cr] <- rep(100, n)
    biomass2[a, t, cr] <- sum(N2[, a, t, cr] * W)
    SSB2[a, t, cr] <- sum(N2[, a, t, cr]*W*Mat)
    abundance_all2[a, t, cr] <- sum(N2[, a, t, cr])
    abundance_mature2[a, t, cr] <- sum(N2[m:n, a, t, cr])
  }

  # Step population forward in time with set fishing level
  for (t in (rec_age + 1):(eq_time - 1)) {
    
    # biology
    PD <- pop_dynamics(a = 1, t, cr = 1, rec_age, max_age, n, SSB2, N2, W, Mat, 
                       A, R0, h, B0, Eps2, sigma_R, Fb, E2, S, M, FM2, m, 
                       abundance_all2, abundance_mature2, biomass2, fishing = F)
    
    SSB2                <- PD[[1]]
    N2                  <- PD[[3]]
    abundance_all2      <- PD[[4]]
    abundance_mature2   <- PD[[5]]
    biomass2            <- PD[[6]]
    
  }
  
  # plotting for troubleshooting
  # plot(1:eq_time, N2[1, 1, 1:eq_time, 1], type = 'l', ylim = c(0, 2e4), col = 'green')
  # for (x in 2:(n - 1)) {
  #   lines(1:eq_time, N2[x, 1, 1:eq_time, 1], col = 'red')
  # }
  # lines(1:eq_time, N2[n, 1, 1:eq_time, 1], col = 'blue')
  
  SAD <- N2[, 1, eq_time - 1, 1]
  
  # # Adjust stable age distribution to give total biomass of B0
  # B1 <- as.numeric(SAD %*% W)
  # SAD1 <- SAD * (B0 / B1)
  
  return(SAD)
  
}