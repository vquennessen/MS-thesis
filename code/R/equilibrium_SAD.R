equilibrium_SAD <- function(rec_age, max_age, n, W, R0,
                            Mat, h, B0, sigma_R, Fb, S, M, eq_time, m, 
                            stochasticity, rho_R, nat_mortality) {
  
  # Initialize population size array
  # Dimensions = age * 1 * time * 1 * 1
  N2 <- array(rep(0, n*eq_time), c(n, 1, eq_time, 1, 1))
  FM2 <- array(rep(0, n*eq_time), c(n, 1, eq_time, 1, 1))
  
  # Initialize biomass and SSB arrays
  # Dimensions = 1 * time * 1 * 1
  biomass2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
  SSB2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
  abundance_all2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
  abundance_mature2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
  
  # Recruitment normal variable
  # Dimensions = area * timeT * CR * 1
  if (stochasticity == T) {
    nuR2 <- array(rnorm(eq_time, 0, sigma_R), c(1, eq_time, 1, 1))
  } else if (stochasticity == F) {
    nuR2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
  }
  
  # Recruitment error
  # Dimensions = area * timeT * CR * 1
  Eps2 <- epsilon(A = 1, eq_time, CR = 1, NM = 1, nuR2, rho_R)
  
  # Start each age class with 10 individuals
  # Enter FM, N, abundance, and biomasses for time = 1 to rec_age
  # Dimensions = age * area * time * CR * M values (3)
  for (t in 1:(rec_age + 1)) {
    N2[, 1, t, 1, 1] <- rep(100, n)
    biomass2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1] * W)
    SSB2[1] <- sum(N2[, 1, t, 1, 1]*W*Mat)
    abundance_all2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1])
    abundance_mature2[1, t, 1, 1] <- sum(N2[m:n, 1, t, 1, 1])
  }
  
  # Step population forward in time with set fishing level
  for (t in (rec_age + 1):(eq_time - 1)) {
    
    # biology
    PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, rec_age, max_age, n, SSB2, 
                       N2, W, Mat, A, R0, h, B0, Eps2, sigma_R, Fb, E2, S, 
                       NM, FM2, m, abundance_all2, abundance_mature2, 
                       biomass2, fishing = F, nat_mortality, recruitment_mode)
    
    FM2[, 1, t, 1, 1]               <- PD[[1]]
    N2[, 1, t, 1, 1]                <- PD[[2]]
    abundance_all2[1, t, 1, 1]      <- PD[[3]]
    abundance_mature2[1, t, 1, 1]   <- PD[[4]]
    biomass2[1, t, 1, 1]            <- PD[[5]]
    SSB2[1, t, 1, 1]                <- PD[[6]]
    
  }
  
  # plotting for troubleshooting
  # plot(1:eq_time, N2[1, 1, 1:eq_time, 1], type = 'l', ylim = c(0, 2e4), col = 'green')
  # for (x in 2:(n - 1)) {
  #   lines(1:eq_time, N2[x, 1, 1:eq_time, 1], col = 'red')
  # }
  # lines(1:eq_time, N2[n, 1, 1:eq_time, 1], col = 'blue')
  
  SAD <- N2[, 1, eq_time - 1, 1, 1]
  
  # # Adjust stable age distribution to give total biomass of B0
  # B1 <- as.numeric(SAD %*% W)
  # SAD1 <- SAD * (B0 / B1)
  
  return(SAD)
  
}