initialize_arrays <- function(A, time1, time2, R0, rec_age, max_age, L1f, L2f, 
                              Kf, a1f, a2f, af, bf, k_mat, Fb, L50, sigma_R, 
                              rho_R, fleets, alpha, beta, start, F_fin, 
                              L50_up, L50_down, cf, switch, full, x, sp, M, CR, 
                              phi, stochasticity, r, D, transects, h, surveys,
                              fishing, error) {
  
  # total amount of timesteps (years)
  timeT <- time1 + time2            
  
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
  
  # Selectivity at age (updated)
  # Dimensions = 1 * age
  S <- selectivity_at_age(fleets, L, max_age, rec_age, alpha, L50_up, L50_down,
                          F_fin, beta, n, cf, age)
  
  # Maturity at age
  # Dimensions = 1 * age
  Mat <- fraction_mature_at_age(n, k_mat, L, L50)
  
  # Cutoff for maturity
  m <- age[min(which(Mat > 0.5))]
  
  # Range of natural mortalities (low, correct, and high)
  nat_mortality <- c(M - error, M, M + error)
  NM <- length(nat_mortality)
  
  # Initialize age-structured population size matrix
  # Dimensions = age * area * time * CR * M values (3)
  N <- array(rep(0, n*A*timeT*CR*NM), c(n, A, timeT, CR, NM))
  
  # Initialize spawning stock biomass array
  # Dimensions = area * time * cr * M values (3)
  SSB <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM))
  
  # Initialize abundance arrays
  # Dimensions = area * time * CR * M values (3)
  abundance_all <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM))
  abundance_mature <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM))
  
  # Initialize biomass array
  # Dimensions = area * time * CR * M values (3)
  biomass <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM)) 
  
  # Recruitment normal variable
  # Dimensions = area * timeT * CR * M values (3)
  if (stochasticity == T) {
    nuR <- array(rnorm(A*timeT*CR*NM, 0, sigma_R), c(A, timeT, CR, NM))
  } else if (stochasticity == F) {
    nuR <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM))
  }
  
  # Recruitment error
  # Dimensions = area * timeT * CR * M values (3)
  Eps <- epsilon(A, timeT, CR, NM, nuR, rho_R)
  
  # Unfished spawning stock biomass
  B0 <- R0/phi  
  
  # Initialize count array
  # Dimensions = area * time * transects * 2 * CR * M values (3)
  Count <- array(rep(0, A*timeT*transects*2*CR*NM), 
                 c(A, timeT, transects, 2, CR, NM))
  
  # Calculate standard deviation of normal variable
  # Based on Babcock & MacCall (2011): Eq. (15)
  sigma_sp <- sqrt(log(1 + (sp/x)^2))
  
  # Sampling normal variable
  # Dimensions = area * timeT * CR * M values (3)
  if (stochasticity == T) {
    nuS <- array(rnorm(A*timeT*CR*NM, 0, sigma_sp), c(A, timeT, CR, NM))
  } else if (stochasticity == F) {
    nuS <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM))
  }
  
  # Calculate delta - constant of proportionality
  # Based on Babcock & MacCall (2011): Eq. (13)
  Delta <- r / D
  
  # Calculate gamma
  # Based on Babcock & MacCall (2011): Eq. (16)
  Gamma <- x / D
  
  # Initialize fishing mortality rate
  # Dimensions = age * area * time * CR * M values (3)
  FM <- array(rep(0, n*A*timeT*CR*NM), c(n, A, timeT, CR, NM))
  
  # Initialize fishing effort in each area
  # Dimensions = area * time * CR * M values (3)
  E <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM)) 
  
  # Initialize catch-at-age matrix
  # Dimensions = age * area * time * CR * M values (3)
  catch <- array(rep(0, n*A*timeT*CR*NM), c(n, A, timeT, CR, NM))
  
  # Initialize yield matrix
  # Dimensions = area * time * CR * M values (3)
  yield <- array(rep(0, A*timeT*CR*NM), c(A, timeT, CR, NM))
  
  if (fishing == T) {
    
    # Initial fishing effort
    E[, 1:time1, , ] <- rep(1/A, A*CR*time1*NM)
    
    # Set constant fishing mortality rate for first 50 years
    fm <- fishing_mortality(a = 1, t = 1, cr = 1, nm = 2, FM, A, Fb, E, S)
    FM[, , 1:time1, , ] <- rep(fm, A*time1*CR*NM)
    
  } 
  
  # Stable age distributions, derived from equilibrium conditions with Fb
  # Dimensions age
  SAD <- equilibrium_SAD(rec_age, max_age, n, W, R0, Mat, h, B0, sigma_R, Fb, 
                         S, M, eq_time = 150, m, stochasticity = F, rho_R, 
                         nat_mortality = M)
  
  # Enter N, abundance, and biomasses for time = 1 to rec_age
  # Dimensions = age * area * time * CR
  for (a in 1:A) {
    for (t in 1:rec_age) {
      for (cr in 1:CR) {
        for (nm in 1:NM) {
        N[, a, t, cr, nm] <- SAD
        abundance_all[a, t, cr, nm] <- sum(N[, a, t, cr, nm])
        abundance_mature[a, t, cr, nm] <- sum(N[m:(max_age-rec_age + 1), a, t, cr, nm])
        biomass[a, t, cr, nm] <- sum(N[, a, t, cr, nm] * W)
        SSB[a, t, cr, nm] <- sum(N[, a, t, cr, nm]*W*Mat)
        }
      }
    }
  }
  
  # initialize relative biomass matrix
  # Dimensions = area * time2 + 1 * CR 
  rel_biomass <- array(rep(0, A*(time2 + 1)*CR), c(A, time2 + 1, CR))
  
  # initialize relative yield matrix
  # Dimensions = area * time2 + 1 * CR 
  rel_yield <- array(rep(0, A*(time2 + 1)*CR), c(A, time2 + 1, CR))
  
  # initialize relative spawning stock biomass matrix
  # Dimensions = area * time2 + 1 * CR
  rel_SSB <- array(rep(0, A*(time2 + 1)*CR), c(A, time2 + 1, CR))
  
  # initialize density ratio matrix
  # Dimensions = timeT * CR
  Density_Ratios <- array(rep(0, (time2 + 1)*CR), c(time2 + 1, CR))
  
  output <- list(timeT, E, n, L, W, Mat, m, S, FM, N, SSB, 
                 abundance_all, abundance_mature, biomass, Count, nuS, 
                 Eps, catch, yield, B0, Delta, Gamma, rel_biomass, rel_yield,
                 rel_SSB, nat_mortality, NM, Density_Ratios)
  
  return(output)
  
}