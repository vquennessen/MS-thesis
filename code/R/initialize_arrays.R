initialize_arrays <- function(A, time, time2, R0, rec_age, max_age, L1f, L2f, 
                              Kf, a1f, a2f, af, bf, k_mat, Fb, L50, sigma_R, 
                              rho_R, fleets, alpha, beta, start, F_fin, 
                              L_50_up, L50_down, cf, switch, full, x, sp, M, CR, 
                              phi, Init_size) {
  
  # total amount of timesteps (years)
  timeT <- time + time2            
  
  # Initialize fishing effort in each area
  # Dimensions = area * time * CR
  E <- array(rep(NaN, A*timeT*CR), c(A, timeT, CR))  
  
  # Initial fishing effort
  E[, 1:time, ] <- rep(1/A, A*CR*time)
  
  # ages for which fish have recruited
  age <- rec_age:max_age 
  
  # number of age classes
  n <- length(age)             
  
  # Length at age
  # Dimensions = 1 * age
  L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f)
  
  # Weight at age
  # Dimensions = 1 * age
  W <- weight_at_age(L, af, bf)
  
  # Maturity at age
  # Dimensions = 1 * age
  Mat <- fraction_mature_at_age(n, k_mat, L, L50)
  
  # Cutoff for maturity
  m <- age[min(which(Mat > 0.5))]
  
  # Selectivity at age
  # Dimensions = 1 * age
  S <- selectivity_at_age(L, fleets, alpha, beta, start, F_fin, L50_up, 
                          L50_down, cf, switch, full)
  
  # Fishing mortality
  # Initialize array
  # Dimensions = age * area * time * CR
  FM <- array(rep(0, n*A*timeT*CR), c(n, A, time + time2, CR))
  
  # Initialize age-structured population size matrix
  # Dimensions = age * area * time * CR
  N <- array(rep(0, n*A*timeT*CR), c(n, A, timeT, CR))
  
  # Initialize spawning stock biomass array
  # Dimensions = area * time
  SSB <- array(rep(0, A*timeT*CR), c(A, timeT, CR))
  
  # Initialize abundance arrays
  # Dimensions = area * time * CR
  abundance_all <- array(rep(0, A*timeT*CR), c(A, timeT, CR))
  abundance_mature <- array(rep(0, A*timeT*CR), c(A, timeT, CR))
  
  # Initialize biomass array
  # Dimensions = area * time * CR
  biomass <- array(rep(0, A*timeT*CR), c(A, timeT, CR)) 
  
  # Initialize count array
  # Dimensions = area * time * transects * 2 * CR
  count_sp <- array(rep(0, A*timeT*transects*2*CR), c(A, timeT, transects, 2, CR))
  
  # Calculate standard deviation of normal variable
  # Based on Babcock & MacCall (2011): Eq. (15)
  sigma_sp <- sqrt(log(1 + (sp/x)^2))
  
  # Sampling normal variable
  # Dimensions = area * timeT * CR
  # nuS <- array(rnorm(A*timeT*CR, 0, sigma_sp), c(A, timeT, CR))
  nuS <- array(rep(0.5, A*timeT*CR), c(A, timeT, CR))

  # Recruitment normal variable
  # nuR <- array(rnorm(A*timeT*CR, 0, sigma_R), c(A, timeT, CR))
  nuR <- array(rep(0.5, A*timeT*CR), c(A, timeT, CR))
  
  # Recruitment error
  # Dimensions = area * timeT * CR
  Eps <- epsilon(A, timeT, CR, nuR, sigma_R, rho_R)
  
  # Length at age for stable age distribution
  # Dimensions = 1 * age (0 to max_age)
  L0 <- length_at_age(0:max_age, L1f, L2f, Kf, a1f, a2f)
  
  # Weight at age for stable age distribution
  # Dimensions = 1 * age (0 to max_age)
  W0 <- weight_at_age(L0, af, bf)
  
  # Stable age distribution
  # Dimensions = 1 * age (0 to max_age)
  SAD <- stable_age_distribution(b, c, max_age, m, L0, W0, rec_age, M, Fb)
  
  # Enter FM, N, abundance, and biomasses for time = 1 to rec_age
  # Dimensions = age * area * time * CR
  for (a in 1:A) {
    for (t in 1:rec_age) {
      for (cr in 1:CR) {
        FM <- fishing_mortality(a, t, cr, FM, A, Fb, E, S)
        N[, a, t, cr] <- Init_size*SAD[rec_age:max_age]
        abundance_all[a, t, cr] <- sum(N[, a, t, cr])
        abundance_mature[a, t, cr] <- sum(N[m:(max_age-1), a, t, cr])
        biomass[a, t, cr] <- sum(N[, a, t, cr] * W)
      }
    }
  }
  
  # Initialize catch-at-age matrix
  # Dimensions = age * area * time * CR
  catch <- array(rep(0, n*A*timeT*CR), c(n, A, timeT, CR))
  
  # Initialize yield matrix
  # Dimensions = area * time * CR
  yield <- array(rep(0, A*timeT*CR), c(A, timeT, CR))
  
  # Unfished spawning stock biomass
  B0 <- R0/phi
  
  output <- list(timeT, E, age, n, L, W, Mat, m, S, FM, N, SSB, 
                 abundance_all, abundance_mature, biomass, count_sp, nuS, 
                 Eps, L0, W0, catch, yield, B0)
  
  return(output)
  
}