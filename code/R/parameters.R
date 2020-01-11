parameters = function(Species) {
  
  # Error handling
  
  if (!is.character(Species)) {
    stop("Study species must be a character string.")
  }
  
  if (Species == 'BR2003') {
    # Black Rockfish 2003 assessment
    # Source: Ralston & Dick 2003
    Max_age <- 35                           # maximum age
    M <- 0.14                               # natural mortality
    Rec_age <- 2                            # age at recruitment
    A <- 1.68e-5; B <- 3                    # weight at length parameters (f)
    A1 <- 5; L1 <- 32.21;                   # growth parameters (f)
    A2 <- 15; L2 <- 47.95; K <- 0.2022  
    L50 <- 39.53                            # length at 50% maturity
    K_mat <- -0.4103                        # slope of maturity curve
    LDP <- 0.1                              # larval drift proportion
    H <- 0.65                               # steepness
    Phi <- 1.1                              # unfished recruits per spawner
    Sigma_R <- 0.5                          # recruitment standard deviation
    Rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.488                              # depletion
    Fb <- 0.2                               # fishing mortality to cause D
    R <- 0.77                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    X <- 15.42                              # mean of positive transects
    SP <- 16.97                             # std of positive transects

    C <- 0.2747                             # eggs produced per kg, intercept
    D <- 0.0941                             # eggs produced per kg, slope

                                            #### selectivity parameters ####
    Fleets <- c('sport', 'hook', 'trawl')   # names of fleets
    Alpha <- c(0.33, 0.6, 0.64)             # slope of upcurve per fleet
    Beta <- c(1.2, 0.6, 0)                  # slope of downcurve per fleet
    F_fin <- c(0.25, 0.06, 1)               # final selectivity if dome-shaped 
    A50_up <- c(2, 5, 10)                   # A50 value for upcurve
    A50_down <- c(6, 16, 35)                # A50 value for downcurve
    Cf <- c(0.71, 0.28, 0.01)               # fraction of fishery caught / fleet
                                            #       from upcurve to 1
                                            
  }
  
  if (Species == 'BR2015') {
    # source: Cope et al. 2016
    Max_age <- 40                           # maximum age
    M <- 0.17                               # natural mortality
    Rec_age <- 3                            # age at recruitment
    A <- 2.6e-5; B <- 2.88                  # weight at length parameters (f)
    A1 <- 1; L1 <- 20.32;                   # growth parameters (f)
    A2 <- 40; L2 <- 49.67; K <- 0.21  
    L50 <- 43.69                            # length at 50% maturity
    K_mat <- -0.66                          # slope of maturity curve
    LDP <- 0.1                              # larval drift proportion
    H <- 0.77                               # steepness
    Phi <- 1.1                              # unfished recruits per spawner
    Sigma_R <- 0.5                          # recruitment standard deviation
    Rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.604                              # depletion
    Fb <- 0.02                              # fishing mortality to cause D 
    R <- 0.77                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    X <- 15.42                              # mean of positive transects
    SP <- 16.97                             # std of positive transects
    
    C <- 0.2747                             # eggs produced per kg, intercept
    D <- 0.0941                             # eggs produced per kg, slope
    
    #### selectivity parameters ####
    Fleets <- c('trawl', 'live', 'dead',    # names of fleets
                'ocean', 'shore')   
    Alpha <- c(0.325, 0.4, 0.35, 
               0.65, 0.425)                 # slope of upcurve per fleet
    Beta <- c(0.25, 0.5, 0.4, 1.1, 0.5)     # slope of downcurve per fleet
    F_fin <- c(0.325, 0.05, -0.11, 
               -0.025, 0.135)               # final selectivity if dome-shaped 
    A50_up <- c(7, 5, 5, 5, 3)              # L50 value for upcurve
    A50_down <- c(15, 13, 13, 12, 6.5)      # L50 value for downcurve
    Cf <- c(0.0001, 0.1679, 0.0982,         # fraction of fishery caught / fleet
            0.6979, 0.0358)               
  } 
  
  if (Species == 'CAB2005') {
    # Black Rockfish 2003 assessment
    # Source: Ralston & Dick 2003
    Max_age <- 17                           # maximum age
    M <- 0.275                              # natural mortality
    Rec_age <- 1                            # age at recruitment
    A <- 1.24e-5; B <- 3.113                # weight at length parameters (f)
    A1 <- 5; L1 <- 41.3;                    # growth parameters (f)
    A2 <- 30; L2 <- 61.9; K <- 0.18  
    L50 <- 25.702                           # length at 50% maturity
    K_mat <- -0.743                         # slope of maturity curve
    LDP <- 0.1                              # larval drift proportion
    H <- 0.7                                # steepness
    Phi <- 0.71                             # unfished recruits per spawner
    Sigma_R <- 1                            # recruitment standard deviation
    Rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.4                                # depletion
    Fb <- 1.1                               # fishing mortality to cause D
    R <- 0.27                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    X <- 3.27                               # mean of positive transects
    SP <- 3.32                              # std of positive transects
    
    C <- 0.0273                             # eggs produced per kg, intercept
    D <- 1.53e-4                            # eggs produced per kg, slope
    
    #### selectivity parameters ####
    Fleets <- c('dead', 'live', 'man.made', 
                'shore', 'PBR', 'CPFV')     # names of fleets
    Alpha <- c(0.4, 0.25, 5, 0.15, 
               0.35, 0.25)                  # slope of upcurve per fleet
    Beta <- c(0, 0.35, 0.5, 0.26, 0, 0)     # slope of downcurve per fleet
    F_fin <- c(1, 0.7, 0.38, 0.43, 1, 1)    # final selectivity if dome-shaped 
    A50_up <- c(2, 3, 1, 1, 2, 1)           # A50 value for upcurve
    A50_down <- c(1, 17, 5, 4, 1, 1)        # A50 value for downcurve
    Cf <- c(0.0186, 0.3389, 0.1208, 0.0241, 
            0.1971, 0.3006)                 # fraction of fishery caught / fleet

  }
  
  if (Species == 'CAB2019') {
    # Black Rockfish 2003 assessment
    # Source: Ralston & Dick 2003
    Max_age <- 20                           # maximum age
    M <- 0.26                               # natural mortality
    Rec_age <- 4                            # age at recruitment
    A <- 1.90e-5; B <- 2.99                 # weight at length parameters (f)
    A1 <- -1.28; L1 <- -0.37;               # growth parameters (f)
    A2 <- 3.86; L2 <- 43.67; K <- 0.225  
    L50 <- 43                               # length at 50% maturity
    K_mat <- -0.7                           # slope of maturity curve
    LDP <- 0.1                              # larval drift proportion
    H <- 0.7                                # steepness
    Phi <- 1.1                              # unfished recruits per spawner
    Sigma_R <- 0.5                          # recruitment standard deviation
    Rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.528                              # depletion
    Fb <-                                # fishing mortality to cause D
    R <- 0.247                              # Proportion of positive transects 
                                            #       in PISCO monitoring data
    X <-                                # mean of positive transects
    SP <-                               # std of positive transects
    
    C <- 0.0273                             # eggs produced per kg, intercept
    D <- 1.53e-4                            # eggs produced per kg, slope
    
    #### selectivity parameters ####
    Fleets <- c('live', 'dead', 'ocean', 'shore') # names of fleets
    Alpha <- c(0.4, 0.33, 0.35, 0.9)              # slope of upcurve per fleet
    Beta <- c(0.35, 0, 0, 0.2)                    # slope of downcurve per fleet
    F_fin <- c(0.7, 1, 1, 0.07)                   # final select. if dome-shaped 
    A50_up <- c(3, 4, 2, 1)                       # A50 value for upcurve
    A50_down <- c(17, 1, 1, 3)                    # A50 value for downcurve
    Cf <- c(0.6033, 0.0415, 0.3423, 0.0130)       # fraction of fishery
    
  }
  
  if (Species == 'LING2017') {
    # Black Rockfish 2003 assessment
    # Source: Ralston & Dick 2003
    Max_age <- 25                           # maximum age
    M <- 0.28                               # natural mortality
    Rec_age <- 3                            # age at recruitment
    A <- 2.76e-6; B <- 3.28                 # weight at length parameters (f)
    A1 <- 1; L1 <- 17.28;                   # growth parameters (f)
    A2 <- 20; L2 <- 120; K <- 0.128  
    L50 <- 56.7                             # length at 50% maturity
    K_mat <- -0.27                          # slope of maturity curve
    LDP <- 0.1                              # larval drift proportion
    H <- 0.7                                # steepness
    Phi <- 1.1                              # unfished recruits per spawner
    Sigma_R <- 0.55                         # recruitment standard deviation
    Rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.579                              # depletion
    Fb <-                                # fishing mortality to cause D
    R <-                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    X <-                                # mean of positive transects
    SP <-                               # std of positive transects
      
    C <- -7.46                              # eggs produced per kg, intercept
    D <- 1.335                              # eggs produced per kg, slope
    
    #### selectivity parameters ####
    Fleets <- c('trawl', 'fixed_gear', 'WArec', 'ORrec') # names of fleets
    Alpha <- c(0.25, 0.25, 0.55, 1)         # slope of upcurve per fleet
    Beta <- c(0.09, 0.3, 0.17, 0.15)        # slope of downcurve per fleet
    F_fin <- c(0.07, 0, 0, 0)               # final select. if dome-shaped 
    A50_up <- c(3, 5, 5, 3)                 # A50 value for upcurve
    A50_down <- c(15, 12, 10, 9)            # A50 value for downcurve
    Cf <- c(0.2872, 0.1379, 0.3253, 0.2496) # fraction of fishery
    
  }
  
  output = list(Max_age, M, Rec_age, A, B, A1, L1f, A2, L2, K, L50, K_mat, LDP, 
                H, Phi, Sigma_R, Rho_R, AMP, D, Fb, R, X, SP, C, D, Fleets, 
                Alpha, Beta, F_fin, A50_up, A50_down, Cf)
  
  return(output)
  
}