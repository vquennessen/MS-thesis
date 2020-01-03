#' Loads life history characteristic parameters for selected species
#'
#' @param species a character string
#' 
#' @return a list of species life history characteristic parameters


parameters = function(species) {
  
  # Error handling
  
  if (!is.character(species)) {
    stop("Study species must be a character string.")
  }
  
  if (species == 'BR2003') {
    # Black Rockfish 2003 assessment
    # Source: Ralston & Dick 2003
    max_age <- 35                           # maximum age
    M <- 0.14                               # natural mortality
    rec_age <- 2                            # age at recruitment
    af <- 1.68e-5; bf <- 3                  # weight at length parameters (f)
    a1f <- 5; L1f <- 32.21;                 # growth parameters (f)
    a2f <- 15; L2f <- 47.95; Kf <- 0.2022  
    L50 <- 39.53                            # length at 50% maturity
    k_mat <- -0.4103                        # slope of maturity curve
    ldp <- 0.1                              # larval drift proportion
    h <- 0.65                               # steepness
    phi <- 1.1                              # unfished recruits per spawner
    sigma_R <- 0.5                          # recruitment standard deviation
    rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.488                              # depletion
    Fb <- 0.2                               # fishing mortality to cause D
    r <- 0.77                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    x <- 15.42                              # mean of positive transects
    sp <- 16.97                             # std of positive transects

    c <- 0.2747                             # eggs produced per kg, intercept
    b <- 0.0941                             # eggs produced per kg, slope

                                            #### fleets: sport, hook, trawl ####
    fleets <- c('sport', 'hook', 'trawl')   # names of fleets
    alpha <- c(0.33, 0.6, 0.64)             # slope of upcurve per fleet
    beta <- c(1.2, 0.6, 0)                  # slope of downcurve per fleet
    F_fin <- c(0.25, 0.06, 0)               # final selectivity if dome-shaped 
    A50_up <- c(2, 5, 10)                   # A50 value for upcurve
    A50_down <- c(6, 16, 35)                # A50 value for downcurve
    cf <- c(0.71, 0.28, 0.01)               # fraction of fishery caught / fleet
                                            #       from upcurve to 1
                                            
  }
  
  if (species == 'BR2015') {
    # source: Cope et al. 2016
    max_age <- 40                           # maximum age
    M <- 0.17                               # natural mortality
    rec_age <- 3                            # age at recruitment
    af <- 2.6e-5; bf <- 2.88                # weight at length parameters (f)
    a1f <- 1; L1f <- 20.32;                 # growth parameters (f)
    a2f <- 40; L2f <- 49.67; Kf <- 0.21  
    L50 <- 43.69                            # length at 50% maturity
    k_mat <- -0.66                          # slope of maturity curve
    ldp <- 0.1                              # larval drift proportion
    h <- 0.77                               # steepness
    phi <- 1.1                              # unfished recruits per spawner
    sigma_R <- 0.5                          # recruitment standard deviation
    rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.604                              # depletion
    Fb <- 0.02                            # fishing mortality to cause D 
                                         ### TODO: find true value of Fb???
    r <- 0.77                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    x <- 15.42                              # mean of positive transects
    sp <- 16.97                             # std of positive transects
    
    c <- 0.2747                             # eggs produced per kg, intercept
    b <- 0.0941                             # eggs produced per kg, slope
    
    #### fleets: sport, hook, trawl ####
    fleets <- c('trawl', 'live', 'dead',    # names of fleets
                'ocean', 'shore')   
    alpha <- c(0.325, 0.4, 0.35, 
               0.65, 0.425)                 # slope of upcurve per fleet
    beta <- c(0.25, 0.5, 0.4, 1.1, 0.5)     # slope of downcurve per fleet
    F_fin <- c(0.325, 0.05, -0.11, 
               -0.025, 0.135)               # final selectivity if dome-shaped 
    A50_up <- c(7, 5, 5, 5, 3)              # L50 value for upcurve
    A50_down <- c(15, 13, 13, 12, 6.5)      # L50 value for downcurve
    cf <- c(0.0001, 0.1679, 0.0982,         # fraction of fishery caught / fleet
            0.6979, 0.0358)               
  } 
  
  if (species == 'CAB2005') {
    # Black Rockfish 2003 assessment
    # Source: Ralston & Dick 2003
    max_age <- 17                           # maximum age
    M <- 0.275                              # natural mortality
    rec_age <- 1                            # age at recruitment
    af <- 1.24e-5; bf <- 3.113              # weight at length parameters (f)
    a1f <- 5; L1f <- 41.3;                  # growth parameters (f)
    a2f <- 30; L2f <- 61.9; Kf <- 0.18  
    L50 <- 25.702                           # length at 50% maturity
    k_mat <- -0.743                         # slope of maturity curve
    ldp <- 0.1                              # larval drift proportion
    h <- 0.7                                # steepness
    phi <- 0.71                             # unfished recruits per spawner
    sigma_R <- 1                            # recruitment standard deviation
    rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.4                                # depletion
    Fb <- 1.1                               # fishing mortality to cause D
    r <- 0.27                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    x <- 3.27                               # mean of positive transects
    sp <- 3.32                              # std of positive transects
    
    c <- 0.0273                             # eggs produced per kg, intercept
    b <- 1.53e-4                            # eggs produced per kg, slope
    
    #### fleets: sport, hook, trawl ####
    fleets <- c('dead', 'live', 'man.made', 
                'shore', 'PBR', 'CPFV')     # names of fleets
    alpha <- c(0.4, 0.25, 5, 0.15, 
               0.35, 0.25)                  # slope of upcurve per fleet
    beta <- c(0, 0.35, 0.5, 0.26, 0, 0)     # slope of downcurve per fleet
    F_fin <- c(1, 0.7, 0.38, 0.43, 1, 1)    # final selectivity if dome-shaped 
    A50_up <- c(2, 3, 1, 1, 2, 1)           # A50 value for upcurve
    A50_down <- c(1, 17, 5, 4, 1, 1)        # A50 value for downcurve
    cf <- c(0.0186, 0.3389, 0.1208, 0.0241, 
            0.1971, 0.3006)                 # fraction of fishery caught / fleet

  }
  
  if (species == 'CAB2019') {
    # Black Rockfish 2003 assessment
    # Source: Ralston & Dick 2003
    max_age <- 20                           # maximum age
    M <- 0.26                               # natural mortality
    rec_age <- 4                            # age at recruitment
    af <- 1.90e-5; bf <- 2.99               # weight at length parameters (f)
    a1f <- -1.28; L1f <- -0.37;             # growth parameters (f)
    a2f <- 3.86; L2f <- 43.67; Kf <- 0.225  
    L50 <- 43                               # length at 50% maturity
    k_mat <- -0.7                           # slope of maturity curve
    ldp <- 0.1                              # larval drift proportion
    h <- 0.7                                # steepness
    phi <- 1.1                              # unfished recruits per spawner
    sigma_R <- 0.5                          # recruitment standard deviation
    rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.528                              # depletion
    Fb <-                                # fishing mortality to cause D
    r <- 0.247                              # Proportion of positive transects 
                                            #       in PISCO monitoring data
    x <-                                # mean of positive transects
    sp <-                               # std of positive transects
    
    c <- 0.0273                             # eggs produced per kg, intercept
    b <- 1.53e-4                            # eggs produced per kg, slope
    
    #### fleets: sport, hook, trawl ####
    fleets <- c('dead', 'live', 'man.made', 
                'shore', 'PBR', 'CPFV')     # names of fleets
    alpha <- c()                  # slope of upcurve per fleet
    beta <- c()     # slope of downcurve per fleet
    F_fin <- c()    # final selectivity if dome-shaped 
    A50_up <- c()           # A50 value for upcurve
    A50_down <- c()        # A50 value for downcurve
    cf <- c()                 # fraction of fishery caught / fleet
    
  }
  
  
  output = list(max_age, M, rec_age, af, bf, a1f, L1f, a2f, L2f, Kf, L50, k_mat, 
                ldp, h, phi, sigma_R, rho_R, AMP, D, Fb, r, x, sp, c, b, fleets, 
                alpha, beta, F_fin, A50_up, A50_down, cf)
  
  return(output)
  
}