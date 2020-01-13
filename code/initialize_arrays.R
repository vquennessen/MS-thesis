initialize_arrays <- function(A, Time1, Time2, R0, Rec_age, Max_age, L1, L2, K, 
                              A1, A2, WA, WB, K_mat, Fb, L50, Sigma_R, Rho_R, 
                              Fleets, Alpha, Beta, Start, F_fin, L50_up, 
                              L50_down, Cf, X, SP, M, CR, Phi, Stochasticity, R, 
                              D, Transects, H, Surveys, Fishing, Error, 
                              Recruitment_mode) {
  
  # total amount of timesteps (years)
  TimeT <- Time1 + Time2            
  
  # ages for which fish have recruited
  Age <- Rec_age:Max_age 
  
  # number of age classes
  Num <- length(Age)             
  
  # Length at age
  # Dimensions = 1 * age
  L <- length_at_age(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = F)
  
  # Weight at age
  # Dimensions = 1 * age
  W <- weight_at_age(L, WA, WB)
  
  # Selectivity at age (updated)
  # Dimensions = 1 * age
  S <- selectivity_at_age(Fleets, L, Max_age, Rec_age, Alpha, A50_up, 
                          A50_down, F_fin, Beta, Num, Cf, Age)
  
  # Maturity at age
  # Dimensions = 1 * age
  Mat <- fraction_mature_at_age(Num, K_mat, L, L50)
  
  # Cutoff for maturity
  A50_mat <- age[min(which(Mat > 0.5))]
  
  # Range of natural mortalities (low, correct, and high)
  Nat_mortality <- c(M - Error, M, M + Error)
  NM <- length(Nat_mortality)
  
  # Initialize age-structured population size matrix
  # Dimensions = age * area * time * CR * M values (3)
  N <- array(rep(0, Num*A*TimeT*CR*NM), c(Num, A, TimeT, CR, NM))
  
  # Initialize spawning stock biomass array
  # Dimensions = area * time * cr * M values (3)
  SSB <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM))
  
  # Initialize abundance arrays
  # Dimensions = area * time * CR * M values (3)
  Abundance_all <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM))
  Abundance_mature <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM))
  
  # Initialize biomass array
  # Dimensions = area * time * CR * M values (3)
  Biomass <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM)) 
  
  # Recruitment normal variable
  # Dimensions = area * timeT * CR * M values (3)
  if (Stochasticity == T) {
    NuR <- array(rnorm(A*TimeT*CR*NM, 0, Sigma_R), c(A, TimeT, CR, NM))
  } else if (Stochasticity == F) {
    NuR <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM))
  }
  
  # Recruitment error
  # Dimensions = area * timeT * CR * M values (3)
  Eps <- epsilon(A, TimeT, CR, NM, NuR, Rho_R)
  
  # Unfished spawning stock biomass
  B0 <- R0/Phi  
  
  # Initialize count array
  # Dimensions = area * time * transects * 2 * CR * M values (3)
  Count <- array(rep(0, A*TimeT*Transects*2*CR*NM), 
                 c(A, TimeT, Transects, 2, CR, NM))
  
  # Calculate standard deviation of normal variable
  # Based on Babcock & MacCall (2011): Eq. (15)
  Sigma_sp <- sqrt(log(1 + (SP/X)^2))
  
  # Sampling normal variable
  # Dimensions = area * timeT * CR * M values (3)
  if (Stochasticity == T) {
    NuS <- array(rnorm(A*TimeT*CR*NM, 0, Sigma_sp), c(A, TimeT, CR, NM))
  } else if (Stochasticity == F) {
    NuS <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM))
  }
  
  # Calculate delta - constant of proportionality
  # Based on Babcock & MacCall (2011): Eq. (13)
  Delta <- R / D
  
  # Calculate gamma
  # Based on Babcock & MacCall (2011): Eq. (16)
  Gamma <- X / D
  
  # Initialize fishing mortality rate
  # Dimensions = age * area * time * CR * M values (3)
  FM <- array(rep(0, Num*A*TimeT*CR*NM), c(Num, A, TimeT, CR, NM))
  
  # Initialize fishing effort in each area
  # Dimensions = area * time * CR * M values (3)
  E <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM)) 
  
  # Initialize catch-at-age matrix
  # Dimensions = age * area * time * CR * M values (3)
  Catch <- array(rep(0, Num*A*TimeT*CR*NM), c(Num, A, TimeT, CR, NM))
  
  # Initialize yield matrix
  # Dimensions = area * time * CR * M values (3)
  Yield <- array(rep(0, A*TimeT*CR*NM), c(A, TimeT, CR, NM))
  
  if (Fishing == T) {
    
    # Initial fishing effort
    E[, 1:Time1, , ] <- rep(1/A, A*CR*Time1*NM)
    
    # Set constant fishing mortality rate for first 50 years
    fm <- fishing_mortality(a = 1, t = 1, cr = 1, nm = 2, FM, A, Fb, E, S)
    FM[, , 1:Time1, , ] <- rep(fm, A*Time1*CR*NM)
    
  } 
  
  # Stable age distributions, derived from equilibrium conditions with Fb
  # Dimensions age
  SAD <- equilibrium_SAD(Rec_age, Max_age, Num, W, R0, Mat, H, B0, Sigma_R, Fb, 
                         S, M, eq_time = 150, A50_mat, Stochasticity = F, Rho_R, 
                         Nat_mortality = M, Recruitment_mode, A)
  
  # Enter N, abundance, and biomasses for time = 1 to rec_age
  # Dimensions = age * area * time * CR
  for (a in 1:A) {
    for (t in 1:Rec_age) {
      for (cr in 1:CR) {
        for (nm in 1:NM) {
        N[, a, t, cr, nm] <- SAD
        Abundance_all[a, t, cr, nm] <- sum(N[, a, t, cr, nm])
        Abundance_mature[a, t, cr, nm] <- sum(N[A50_mat:(Max_age-Rec_age + 1), a, t, cr, nm])
        Biomass[a, t, cr, nm] <- sum(N[, a, t, cr, nm] * W)
        SSB[a, t, cr, nm] <- sum(N[, a, t, cr, nm]*W*Mat)
        }
      }
    }
  }
  
  # initialize relative biomass matrix
  # Dimensions = area * time2 + 1 * CR 
  Rel_biomass <- array(rep(0, A*(Time2 + 1)*CR), c(A, Time2 + 1, CR))
  
  # initialize relative yield matrix
  # Dimensions = area * time2 + 1 * CR 
  Rel_yield <- array(rep(0, A*(Time2 + 1)*CR), c(A, Time2 + 1, CR))
  
  # initialize relative spawning stock biomass matrix
  # Dimensions = area * time2 + 1 * CR
  Rel_SSB <- array(rep(0, A*(Time2 + 1)*CR), c(A, Time2 + 1, CR))
  
  # initialize density ratio matrix
  # Dimensions = timeT * CR
  Density_ratio <- array(rep(0, (Time2 + 1)*CR), c(Time2 + 1, CR))
  
  output <- list(TimeT, E, Num, L, W, Mat, A50_mat, S, FM, N, SSB, 
                 Abundance_all, Abundance_mature, Biomass, Count, NuS, 
                 Eps, Catch, Yield, B0, Delta, Gamma, Rel_biomass, Rel_yield,
                 Rel_SSB, Nat_mortality, NM, Density_ratio, Age)
  
  return(output)
  
}