equilibrium_SAD <- function(a, cr, A, rec_age, max_age, n, W, R0,
                            Mat, h, B0, Eps, sigma_R, Fb, S, M, season,
                            catch_form, eq_time, m, stochasticity, rho_R) {
  
  # Initialize population size and catch arrays
  # Dimensions = age * 1 * time * 1
  N2 <- catch2 <- FM2 <- array(rep(0, n*eq_time), c(n, 1, eq_time, 1))
  
  # Initialize biomass and SSB arrays
  # Dimensions = 1 * time * 1
  biomass2 <- yield2 <- SSB2 <- E2 <- array(rep(0, eq_time), c(1, eq_time, 1))
  abundance_all2 <- abundance_mature2 <- array(rep(0, eq_time), c(1, eq_time, 1))
  
  # Recruitment normal variable
  # Dimensions = area * timeT * CR
  if (stochasticity == T) {
    nuR2 <- array(rnorm(1*eq_time*1, 0, sigma_R), c(1, eq_time, 1))
  } else if (stochasticity == F) {
    nuR2 <- array(rep(0.5, 1*eq_time*1), c(1, eq_time, 1))
  }
  
  # Recruitment error
  # Dimensions = area * timeT * CR
  Eps2 <- epsilon(1, eq_time, 1, nuR2, rho_R)
  
  # Initial fishing effort
  E2[, 1:eq_time, ] <- rep(1/A, eq_time)
  
  # Start each age class with 10 individuals
  # Enter FM, N, abundance, and biomasses for time = 1 to rec_age
  # Dimensions = age * area * time * CR
  for (t in 1:(rec_age + 1)) {
    FM2[, 1, t, 1] <- fishing_mortality(1, t, 1, FM2, A, Fb, E2, S)
    N2[, 1, t, 1] <- rep(10, n)
    biomass2[1, t, 1] <- sum(N2[, 1, t, 1] * W)
    SSB2[1, t, 1] <- spawning_stock_biomass(1, t, 1, rec_age, N2, W, Mat)
    catch2[, 1, t, 1] <- catch_at_age(1, t, 1, FM2, M, N2, A, Fb, E2, catch2, 
                                      catch_form, season)
    yield2[1, t, 1] <- sum(catch2[, a, t, cr]*W)
  }

  # Step population forward in time with set fishing level
  for (t in (rec_age + 1):(eq_time - 1)) {
    
    # biology
    PD <- pop_dynamics(a = 1, t, cr = 1, rec_age, max_age, n, SSB2, N2, W, Mat, 
                       A, R0, h, B0, Eps2, sigma_R, Fb, E2, S, M, FM2, m, 
                       abundance_all2, abundance_mature2, biomass2)
    
    SSB2                <- PD[[1]]
    FM2                 <- PD[[2]]
    N2                  <- PD[[3]]
    biomass2            <- PD[[6]]
    
    # fishing
    catch2[, 1, t, 1] <- catch_at_age(1, t, 1, FM2, M, N2, A, Fb, E2, catch2, 
                                     catch_form, season)
    N2[, 1, t, 1] <- N2[, 1, t, 1] - catch2[, 1, t, 1]
    yield2[1, t, 1] <- sum(catch2[, 1, t, 1]*W)
    
  }
  
  SAD <- N2[, 1, eq_time - 1, 1]
  
  # Adjust stable age distribution to give total biomass of B0
  B1 <- as.numeric(SAD %*% W)
  SAD <- SAD * (B0 / B1)
  
  return(SAD)
  
}