initialize_arrays <- function(A, time1, time2, R0, rec_age, max_age, L1f, L2f, 
                              Kf, a1f, a2f, af, bf, k_mat, Fb, L50, sigma_R, 
                              rho_R, fleets, alpha, beta, start, F_fin, 
                              L_50_up, L50_down, cf, switch, full, x, sp, M, CR, 
                              phi, catch_form, season, stochasticity, r, D) {
  
  # total amount of timesteps (years)
  timeT <- time1 + time2            
  
  # Initialize fishing effort in each area
  # Dimensions = area * time * CR
  E <- array(rep(NaN, A*timeT*CR), c(A, timeT, CR))  
  
  # Initial fishing effort
  E[, 1:time1, ] <- rep(1/A, A*CR*time1)
  
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
  FM <- array(rep(0, n*A*timeT*CR), c(n, A, timeT, CR))
  
  # Initialize age-structured population size matrix
  # Dimensions = age * area * time * CR
  N <- array(rep(0, n*A*timeT*CR), c(n, A, timeT, CR))
  
  # Initialize spawning stock biomass array
  # Dimensions = area * time * cr
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
  if (stochasticity == T) {
    nuS <- array(rnorm(A*timeT*CR, 0, sigma_sp), c(A, timeT, CR))
  } else if (stochasticity == F) {
    nuS <- array(rep(0.5, A*timeT*CR), c(A, timeT, CR))
  }

  # Recruitment normal variable
  # Dimensions = area * timeT * CR
  if (stochasticity == T) {
    nuR <- array(rnorm(A*timeT*CR, 0, sigma_R), c(A, timeT, CR))
  } else if (stochasticity == F) {
    nuR <- array(rep(0.5, A*timeT*CR), c(A, timeT, CR))
  }
  
  # Recruitment error
  # Dimensions = area * timeT * CR
  Eps <- epsilon(A, timeT, CR, nuR, sigma_R, rho_R)
  
  # Unfished spawning stock biomass
  B0 <- R0/phi  
  
  # Initialize catch-at-age matrix
  # Dimensions = age * area * time * CR
  catch <- array(rep(0, n*A*timeT*CR), c(n, A, timeT, CR))
  
  # Initialize yield matrix
  # Dimensions = area * time * CR
  yield <- array(rep(0, A*timeT*CR), c(A, timeT, CR))
  
  # # Length at age for stable age distribution
  # # Dimensions = 1 * age (0 to max_age)
  # L0 <- length_at_age(0:max_age, L1f, L2f, Kf, a1f, a2f)
  # 
  # # Weight at age for stable age distribution
  # # Dimensions = 1 * age (0 to max_age)
  # W0 <- weight_at_age(L0, af, bf)
  # 
  # # Stable age distribution, derived from Leslie matrix
  # # Dimensions = 1 * age (0 to max_age)
  # SAD <- Leslie_SAD(b, c, max_age, m, L0, W0, rec_age, M, Fb, h, R0, W)

  # Stable age distribution, derived from equilibrium conditions with Fb
  eq_time <- 150
  SAD <- equilibrium_SAD(1, 1, allocation, A, rec_age, max_age, n, W, R0,
                              Mat, h, B0, Eps, sigma_R, Fb, S, M, season,
                              catch_form, eq_time, m, stochasticity = F)
  
  # # Initial size of whole population at time = 1, 2
  # Init_size <- initial_size(SAD)

  # Enter FM, N, abundance, and biomasses for time = 1 to rec_age
  # Dimensions = age * area * time * CR
  for (a in 1:A) {
    for (t in 1:rec_age) {
      for (cr in 1:CR) {
        FM[, a, t, cr] <- fishing_mortality(a, t, cr, FM, A, Fb, E, S)
        N[, a, t, cr] <- SAD
        abundance_all[a, t, cr] <- sum(N[, a, t, cr])
        abundance_mature[a, t, cr] <- sum(N[m:(max_age-1), a, t, cr])
        biomass[a, t, cr] <- sum(N[, a, t, cr] * W)
        catch[, a, t, cr] <- catch_at_age(a, t, cr, FM, M, N, A, Fb, E, catch, 
                                          catch_form, season)
        yield[a, t, cr] <- sum(catch[, a, t, cr]*W)
        SSB[a, t, cr] <- spawning_stock_biomass(a, t, cr, N, W, Mat)
      }
    }
  }
  
  # Calculate delta - constant of proportionality
  # Based on Babcock & MacCall (2011): Eq. (13)
  delta <- r / D
  
  output <- list(timeT, E, n, L, W, Mat, m, S, FM, N, SSB, 
                 abundance_all, abundance_mature, biomass, count_sp, nuS, 
                 Eps, catch, yield, B0, delta)
  
  return(output)
  
}