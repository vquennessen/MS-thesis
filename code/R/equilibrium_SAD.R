equilibrium_SAD <- function(a, cr, allocation, A, E, rec_age, max_age, n, W, R0,
                            Mat, h, B0, Eps, sigma_R, Fb, E, S, M, FM, season,
                            catch_form, equilibrium_time, m) {
  
  # Initialize population size and catch arrays
  # Dimensions = age * 1 * time * 1
  N <- catch <- array(rep(0, n*equilibrium_time), c(n, 1, equilibrium_time, 1))
  
  # Initialize biomass and SSB arrays
  # Dimensions = 1 * time * 1
  biomass <- SSB <- array(rep(0, equilibrium_time), c(1, equilibrium_time, 1))
  
  # Start each age class with 10 individuals
  N[, 1, 1, 1] <- rep(10, n)
  biomass[1, 1, 1] <- sum(N[, 1]*W)
  SSB[1, 1, 1] <- sum(N[m:n, 1, 1, 1]*W)
  
  # Step population forward in time with set fishing level
  for (t in 2:equilibrium_time) {
    
    # effort allocation
    E <- effort_allocation(1, t, 1, allocation, A, E, biomass)
    
    # biology
    PD <- pop_dynamics(1, t, 1, rec_age, max_age, n, SSB, N, W, Mat, A, R0, 
                       h, B0, Eps, sigma_R, Fb, E, S, M)
    
    SSB                <- PD[[1]]
    FM                 <- PD[[2]]
    N                  <- PD[[3]]
    biomass            <- PD[[6]]
    
    # fishing
    catch[, 1, t, 1] <- catch_at_age(1, t, 1, FM, M, N, A, Fb, E, catch, 
                                     catch_form, season)
    N[, 1, t, 1] <- N[, 1, t, 1] - catch[, 1, t, 1]
    
  }
  
  SAD <- N[, 1, equilibrium_time, 1]
  
  return(SAD)
  
}