initialize_arrays <- function(time, time2, init_effort, rec_age, max_age, L1f, 
                              L2f, Kf, a1f, a2f, af, bf, k_mat, Fb, L50, 
                              sigma_R, rho_R, fleets, alpha, beta, start, F_fin, 
                              L_50_up, L50_down, cf, switch, full, A, x, sp, 
                              initial, M, CR, R0, phi) {
  
  # total amount of timesteps (years)
  timeT <- time + time2            
  
  # nominal fishing effort in each area
  E <- array(rep(init_effort, A), c(1, A))     
  
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
  S <- selectivity_at_age(L, fleets, alpha, beta, start, F_fin, L50_up, L50_down, 
                          cf, switch, full)
  
  # Fishing mortality
  # Initialize array
  # Dimensions = age * area * time
  FM <- array(rep(0, n*A*timeT), c(n, A, time + time2))
  
  # Initial fishing mortality
  FM[, , 1] <- fishing_mortality(A, Fb, E[1], S)
  
  # Initialize age-structured population size matrix
  # Dimensions = age * area * time
  N <- array(rep(0, n*A*timeT), c(n, A, timeT))
  
  # Initialize spawning stock biomass array
  # Dimensions = area * time
  SSB <- array(rep(0, A*timeT), c(A, timeT))
  
  # Initialize recruitment array, 
  # Dimensions = area * time
  R <- array(rep(0, A*timeT), c(A, timeT))
  
  # Initialize abundance arrays
  # Dimensions = area * time
  abundance_all <- array(rep(0, A*timeT), c(A, timeT))
  abundance_mature <- array(rep(0, A*timeT), c(A, timeT))
  
  # Initialize biomass array
  # Dimensions = area * time
  biomass <- array(rep(0, A*timeT), c(A, timeT)) 
  
  # Initialize count array
  # Dimensions = area * time * transects * 2
  time_cnt <- time2 + 3
  count_sp <- array(rep(0, A*time_cnt*transects*2), c(A, time_cnt, transects, 2))
  
  # Enter abundance and biomasses for time = 1, 2
  for (a in 1:A) {
    for (t in 1:2) {
      abundance_all[a, t] <- sum(N[, a, t])
      abundance_mature[a, t] <- sum(N[m:(max_age-1), a, t])
      biomass[a, t] <- sum(N[, a, t] * W)
    }
  }
  
  # Calculate standard deviation of normal variable
  # Based on Babcock & MacCall (2011): Eq. (15)
  sigma_sp <- sqrt(log(1 + (sp/x)^2))
  
  # Sampling normal variable
  # Dimensions = area * time2 + 3 * CR
  # nuS <- array(rnorm(A*(time2 + 3)*CR, 0, sigma_sp), c(A, time2 + 3, CR))
  nuS <- array(rep(0.5, A*(time2 + 3)*CR), c(A, time2 + 3, CR))

  # Recruitment normal variable
  # nuR <- array(rnorm(A*timeT*CR, 0, sigma_R), c(A, timeT, CR))
  nuR <- array(rep(0.5, A*timeT*CR), c(A, timeT, CR))
  
  # Recruitment error
  # Dimensions = area * timeT * CR
  Eps <- epsilon(A, timeT, CR, nuR, sigma_R, rho_R)
  
  #####  Extra vectors to determine stable age distribution
  # Length at age
  # Dimensions = 1 * age (0 to max_age)
  L0 <- length_at_age(0:max_age, L1f, L2f, Kf, a1f, a2f)
  
  # Weight at age
  # Dimensions = 1 * age (0 to max_age)
  W0 <- weight_at_age(L0, af, bf)
  
  # Stable age distribution
  # Dimensions = 1 * age (0 to max_age)
  SAD <- stable_age_distribution(b, c, max_age, m, L0, W0, rec_age, M, Fb)
  
  # Initial age structure
  N[, , 1] <- N[, , 2] <- initial*SAD[rec_age:max_age]
  
  # Enter abundance and biomasses for time = 1, 2
  for (a in 1:A) {
    for (t in 1:2) {
      abundance_all[a, t] <- sum(N[, a, t])
      abundance_mature[a, t] <- sum(N[m:(max_age-1), a, t])
      biomass[a, t] <- sum(N[, a, t] * W)
    }
    
  }
  
  # Initialize catch-at-age matrix
  # Dimensions = age * area * time
  catch <- array(rep(0, n*A*timeT), c(n, A, timeT))
  
  # Initialize yield matrix
  # Dimensions = area * time
  yield <- array(rep(0, A*timeT), c(A, timeT))
  
  # Unfished spawning stock biomass
  B0 <- R0/phi
  
  output <- list(timeT, E, age, n, L, W, Mat, m, S, FM, N, SSB, R, 
                 abundance_all, abundance_mature, biomass, count_sp, nuS, 
                 Eps, L0, W0, catch, yield, B0)
  
  return(output)
  
}