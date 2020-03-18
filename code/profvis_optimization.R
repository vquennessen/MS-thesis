library(profvis)
profvis({
  
  Species = 'BR_OR_2015'
  R0 = 1e5
  A = 5
  MPA = 3
  Time1 = 50
  Time2 = 20
  Error = 0.05
  Transects = 24
  Plotting = FALSE
  Final_DRs = c(0.2, 0.6, 1.0)
  Control_rules = c(1:6)
  Plot_individual_runs <- F
  BM <- FALSE
  Output.Yield = TRUE
  Output.Biomass = TRUE
  Output.SSB = TRUE
  Output.N = TRUE
  Output.Density.Ratio = TRUE
  Output.Effort = TRUE
  Surveys = TRUE
  Fishing = TRUE
  Stochasticity = TRUE
  Recruitment_mode = 'pool'
  Allocation = 'IFD'
  Adult_movement = TRUE
  Fishery_management = TRUE
  Years_sampled = 1
  Areas_sampled = 'all'
  Ind_sampled = 'all'
  
  ##### Load life history characteristics for species ##########################
  
  par <- parameters(Species)
  
  Max_age                <- par[[1]]        # maximum age
  M                      <- par[[2]]        # natural mortality
  Rec_age                <- par[[3]]        # age at recruitment
  WA <- par[[4]];  WB    <- par[[5]]        # weight at length parameters (f)
  A1 <- par[[6]];  L1    <- par[[7]]        # growth parameters (f)
  A2 <- par[[8]];  L2    <- par[[9]]
  K                      <- par[[10]]
  L50                    <- par[[11]]       # length at 50% maturity
  K_mat                  <- par[[12]]       # slope of maturity curve
  LDP                    <- par[[13]]       # larval drift proportion
  H                      <- par[[14]]       # steepness
  Phi                    <- par[[15]]       # unfished recruits per spawner
  Sigma_R                <- par[[16]]       # recruitment standard deviation
  Rho_R                  <- par[[17]]       # recruitment autocorrelation
  AMP                    <- par[[18]]       # adult movement proportion
  D                      <- par[[19]]       # depletion
  Fb                     <- par[[20]]       # fishing mortality to cause D
  P                      <- par[[21]]       # proportion of positive transects
  X                      <- par[[22]]       # mean of positive transects
  SP                     <- par[[23]]       # std of positive transects
  Fleets                 <- par[[24]]       # fishery fleet names
  Alpha                  <- par[[25]]       # slope for upcurve
  Beta                   <- par[[26]]       # slope for downcurve
  F_fin                  <- par[[27]]       # F_fin for fishery, 0 if asymptotic
  A50_up                 <- par[[28]]       # A50 for upcurve
  A50_down               <- par[[29]]       # A50 for downcurve
  Cf                     <- par[[30]]       # fraction of fishery caught / fleet
  
  ##### Population Dynamics - Non-Time Varying #################################
  
  # Initialize arrays for time-varying dynamics
  IA <- initialize_arrays(A, MPA, Final_DRs, Time1, Time2, R0, Rec_age, Max_age,
                          A1, L1, A2, L2, K, WA, WB, K_mat, Fb, L50, Sigma_R,
                          Rho_R, Fleets, Alpha, A50_up, A50_down, F_fin, Beta,
                          Cf, P, X, SP, M, Control_rules, Phi, Stochasticity, D,
                          Transects, H, Surveys, Fishing, Error,
                          Recruitment_mode)
  
  Inside           <- IA[[1]]     # Area(s) in the marine reserve
  Outside          <- IA[[2]]     # Areas not in the marine reserve
  FDR              <- IA[[3]]
  TimeT            <- IA[[4]]     # total amount of timesteps (years)
  L                <- IA[[5]]     # Length at age, dim = 1*age
  W                <- IA[[6]]     # Weight at age, dim = 1*age
  S                <- IA[[7]]     # Selectivity at age
  Mat              <- IA[[8]]     # Fraction mature at age, dim = 1*age
  A50_mat          <- IA[[9]]     # Age at which fraction mature > 0.5
  CR               <- IA[[10]]    # Number of control rules
  Nat_mortality    <- IA[[11]]    # Range of potential natural mortality values
  NM               <- IA[[12]]    # Number of potential natural mortality values
  N                <- IA[[13]]    # Population size, dim = age*area*time
  SSB              <- IA[[14]]    # Spawning stock biomass, dim = area*time
  Abundance_all    <- IA[[15]]    # Abundance, dim = area*time
  Abundance_mature <- IA[[16]]    # Abundance, dim = area*time
  Biomass          <- IA[[17]]    # Biomass, dim = area*time
  Eps              <- IA[[18]]    # Epsilon vector, dim = area*time*CR
  B0               <- IA[[19]]    # Unfished spawning stock biomass
  Count            <- IA[[20]]    # Species count when sampling, dim = area*time
  Sigma_S          <- IA[[21]]    # Sampling normal standard deviation
  NuS              <- IA[[22]]    # Sampling normal variable, dim = area*time*CR
  Delta            <- IA[[23]]    # Constant of proportionality
  Gamma            <- IA[[24]]    # Gamma
  FM               <- IA[[25]]    # Fishing mortality rate, dim = age*area*time
  E                <- IA[[26]]    # nominal fishing effort in each area
  Catch            <- IA[[27]]    # Catch at age
  Yield            <- IA[[28]]    # Yield per area
  Rel_biomass      <- IA[[29]]    # Relative biomass after reserve implementation
  Rel_yield        <- IA[[30]]    # Relative yield after reserve implementation
  Rel_SSB          <- IA[[31]]    # Relative SSB after reserve implementation
  Density_ratio    <- IA[[32]]    # Density ratios
  ENM              <- IA[[33]]    # nm value that represents "true" M
  
  ##### Population Dynamics - Time Varying #####################################
  
  for (fdr in 1:FDR) {
    
    for (t in (Rec_age + 1):(Time1 - 1)) {
      
      for (cr in 1:CR) {
        
        for (nm in 1:NM) {
          
          # effort allocation
          E <- effort_allocation(t, cr, nm, fdr, Allocation, E, Yield, Time1,
                                 Inside, Outside)
          
          # If there is adult movement, add movement
          if (Adult_movement == TRUE) {N <- movement(t, cr, nm, fdr, N, A, AMP)}
          
          # Recruitment / larval movement (if applicable)
          R <- recruitment(t, cr, nm, fdr, SSB, A, R0, H, B0, Eps, Sigma_R,
                           Rec_age, Recruitment_mode, LDP)
          
          for (a in 1:A) {
            
            # biology
            PD <- pop_dynamics(a, t, cr, nm, fdr, Rec_age, Max_age, SSB,
                               N, W, Mat, A, Fb, E, S, NM, FM, A50_mat,
                               Abundance_all, Abundance_mature, Biomass, Fishing,
                               Nat_mortality, R)
            
            FM[, a, t, cr, nm, fdr]               <- PD[[1]]
            N[, a, t, cr, nm, fdr]                <- PD[[2]]
            Abundance_all[a, t, cr, nm, fdr]      <- PD[[3]]
            Abundance_mature[a, t, cr, nm, fdr]   <- PD[[4]]
            Biomass[a, t, cr, nm, fdr]            <- PD[[5]]
            SSB[a, t, cr, nm, fdr]                <- PD[[6]]
            
            # sampling
            if (Surveys == TRUE) {
              Count[a, t, , , cr, nm, fdr] <- sampling(a, t, cr, nm, fdr, Delta,
                                                       Gamma, Abundance_all,
                                                       Abundance_mature,
                                                       Transects, X, Count, NuS,
                                                       A)
            }
            
            # fishing
            if (Fishing == TRUE) {
              Catch[, a, t, cr, nm, fdr] <- catch(a, t, cr, nm, fdr, FM,
                                                  Nat_mortality, N, A, Fb, E,
                                                  Catch)
              Yield[a, t, cr, nm, fdr] <- sum(Catch[, a, t, cr, nm, fdr]*W)
            }
          }
        }
      }
    }
  }
  
  ##### Implement Reserve, and apply control rules #############################
  
  for (fdr in 1:FDR) {
    
    for (t in Time1:TimeT) {
      
      for (cr in 1:CR) {
        
        for (nm in 1:NM) {
          
          # effort allocation
          E <- effort_allocation(t, cr, nm, fdr, Allocation, E, Yield, Time1,
                                 Inside, Outside)
          
          # If there is adult movement, add movement
          if (Adult_movement == TRUE) {N <- movement(t, cr, nm, fdr, N, A, AMP)}
          
          # Recruitment / larval movement (if applicable)
          R <- recruitment(t, cr, nm, fdr, SSB, A, R0, H, B0, Eps, Sigma_R,
                           Rec_age, Recruitment_mode, LDP)
          
          for (a in 1:A) {
            
            # biology
            PD <- pop_dynamics(a, t, cr, nm, fdr, Rec_age, Max_age, SSB,
                               N, W, Mat, A, Fb, E, S, NM, FM, A50_mat,
                               Abundance_all, Abundance_mature, Biomass, Fishing,
                               Nat_mortality, R)
            
            FM[, a, t, cr, nm, fdr]               <- PD[[1]]
            N[, a, t, cr, nm, fdr]                <- PD[[2]]
            Abundance_all[a, t, cr, nm, fdr]      <- PD[[3]]
            Abundance_mature[a, t, cr, nm, fdr]   <- PD[[4]]
            Biomass[a, t, cr, nm, fdr]            <- PD[[5]]
            SSB[a, t, cr, nm, fdr]                <- PD[[6]]
            
            # sampling
            if (Surveys == TRUE) {
              Count[a, t, , , cr, nm, fdr] <- sampling(a, t, cr, nm, fdr, Delta,
                                                       Gamma, Abundance_all,
                                                       Abundance_mature,
                                                       Transects, X, Count, NuS,
                                                       A)
            }
            
            # fishing
            if (Fishing == TRUE) {
              Catch[, a, t, cr, nm, fdr] <- catch(a, t, cr, nm, fdr, FM,
                                                  Nat_mortality, N, A, Fb, E,
                                                  Catch)
              Yield[a, t, cr, nm, fdr] <- sum(Catch[, a, t, cr, nm, fdr]*W)
            }
          }
        }
        
        # management
        if (Fishery_management == TRUE && t < TimeT) {
          E[, t, cr, , fdr] <- control_rule(t, cr, nm, fdr, A, E, Count, Time1,
                                            TimeT, Transects, Nat_mortality,
                                            Final_DRs, Inside, Outside,
                                            Areas_sampled, Ind_sampled,
                                            Years_sampled, BM)
        }
        
        # calculate true density ratio
        Density_ratio <- true_DR(t, cr, fdr, Abundance_all, Inside, Outside,
                                 Density_ratio, Time1, ENM)
        
      }
    }
  }
  
})
