initialize_arrays <- function(L1f, L2f, Kf, a1f, a2f, af, bf, k_mat, Fb,
                              L50, sigma_R, rho_R, fleets, alpha, beta, start, 
                              F_fin, L_50_up, L50_down, cf, switch, full, age, 
                              n, A, time, E, x, sp) {
  
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
  FM <- array(rep(0, n*A*time), c(n, A, time))
  
  # Initial fishing mortality
  FM[, , 1] <- fishing_mortality(A, Fb, E, S)
  
  # Recruitment error
  # Dimensions = 1 * time
  e <- epsilon(time, sigma_R, rho_R)
  
  # Initialize age-structured population size matrix
  # Dimensions = age * area * time
  N <- array(rep(0, n*A*time), c(n, A, time))
  
  # Initial age structure
  initial <- 100          # start with 100 individuals in each area at t = 1
  N[, , 1] <- array(rep(initial, n), c(1, 1, n))
  
  # Initialize spawning stock biomass array
  # Dimensions = area * time
  SSB <- array(rep(0, A*time), c(A, time))
  
  # Initialize recruitment array, 
  # Dimensions = area * time
  R <- array(rep(0, A*time), c(A, time))
  
  # Initialize abundance arrays
  # Dimensions = area * time
  abundance_all <- array(rep(0, A*time), c(A, time))
  abundance_mature <- array(rep(0, A*time), c(A, time))
  
  # Initialize biomass array
  # Dimensions = area * time
  biomass <- array(rep(0, A*time), c(A, time)) 
  
  # Initialize count array
  # Dimensions = area * time * 2
  count_sp <- array(rep(0, A*time*2), c(A, time, 2))
  
  # Calculate standard deviation of normal variable
  # Based on Babcock & MacCall (2011): Eq. (15)
  sigma_sp <- sqrt(log(1 + (sp/x)^2))
  
  # Sampling normal variable
  # Dimensions = area * time
  nu <- array(rnorm(A*time, 0, sigma_sp), c(A, time))
  
  output <- list(L, W, Mat, m, S, FM, e, N, SSB, R, abundance_all, 
                 abundance_mature, biomass, count_sp, nu)
  
  return(output)
  
}