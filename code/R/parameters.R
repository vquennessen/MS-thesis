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
  
  if (species == 'black rockfish 2003') {
    # source: Babcock & MacCall 2011
    max_age <- 35                           # maximum age
    M <- 0.14                               # natural mortality
    rec_age <- 2                            # age at recruitment
    af <- 1.68e-5; bf <- 3                  # weight at length parameters (f)
    am <- 1.68e-5; bm <- 3                  # weight at length parameters (m)
    a1f <- 5; L1f <- 32.21;                 # growth parameters (f)
    a2f <- 15; L2f <- 47.95; Kf <- 0.2022  
    a1m <- 5; L1m <- 31.88;                 # growth parameters (m)
    a2m <- 15; L2m <- 45.39; Km <- 0.1979  
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
    alpha <- c(0.45, 0.369, 0.426)          # slope of upcurve per fleet
    beta <- c(1.28, 0.419, 0)               # slope of downcurve per fleet
    start <- c(20, 20, 20)                  # vulnerability start length / fleet
    F_fin <- c(0.265, 0.29, 0)              # final selectivity if dome-shaped 
    L50_up <- c(25.2, 33.2, 46)             # L50 value for upcurve
    L50_down <- c(36.7, 46.4, 50)           # L50 value for downcurve
    cf <- c(0.71, 0.28, 0.01)               # fraction of fishery caught / fleet
    switch <- c(29, 29, 0)                  # length where selectivity switches 
                                            #       from upcurve to 1
    full <- c(31, 39, 0)                    # length at which downcurve starts
    catch_form <- 'continuous'              # discrete or continuous catch
    season <- 0.5                           # if catch_formulation = discrete, 
                                            #       time at which fishing occurs
                                            #       0 = beginning, 1 = end
                                            
  }
  
  if (species == 'black rockfish 2015') {
    # source: Cope et al. 2016
    max_age <- 40                           # maximum age
    M <- 0.17                               # natural mortality
    rec_age <- 3                            # age at recruitment
    af <- 2.6e-5; bf <- 2.88                # weight at length parameters (f)
    am <- 2.58e-5; bm <- 2.89               # weight at length parameters (m)
    a1f <- 1; L1f <- 20.32;                 # growth parameters (f)
    a2f <- 40; L2f <- 49.67; Kf <- 0.21  
    a1m <- 1; L1m <- 17.47;                 # growth parameters (m)
    a2m <- 40; L2m <- 43.27; Km <- 0.34  
    L50 <- 43.69                            # length at 50% maturity
    k_mat <- -0.66                          # slope of maturity curve
    ldp <- 0.1                              # larval drift proportion
    h <- 0.77                               # steepness
    phi <- 1.1                              # unfished recruits per spawner
    sigma_R <- 0.5                          # recruitment standard deviation
    rho_R <- 0                              # recruitment autocorrelation
    AMP <- 0.1                              # adult movement proportion
    D <- 0.604                              # depletion
    Fb <- 0.08                            # fishing mortality to cause D 
                                         ### TODO: find true value of Fb???
    r <- 0.77                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    x <- 15.42                              # mean of positive transects
    sp <- 16.97                             # std of positive transects
    
    c <- 0.2747                             # eggs produced per kg, intercept
    b <- 0.0941                             # eggs produced per kg, slope
    
    #### fleets: sport, hook, trawl ####
    fleets <- c('trawl', 'live', 'dead',    # names of fleets
                'recO', 'recS')   
    alpha <- c(0.325, 0.4, 0.35, 
               0.65, 0.425)                 # slope of upcurve per fleet
    beta <- c(0.25, 0.5, 0.4, 1.1, 0.5)     # slope of downcurve per fleet
    start <- c(0, 0, 0, 0, 0)
    F_fin <- c(0.325, 0.05, -0.11, 
               -0.025, 0.135)               # final selectivity if dome-shaped 
    L50_up <- c(7, 5, 5, 5, 3)              # L50 value for upcurve
    L50_down <- c(15, 13, 13, 12, 6.5)      # L50 value for downcurve
    cf <- c(0.0001, 0.1679, 0.0982,         # fraction of fishery caught / fleet
            0.6979, 0.0358)               
    switch <- c(0, 0, 0, 0, 0)
    full <- c(0, 0, 0, 0, 0)
    catch_form <- 'continuous'              # discrete or continuous catch
    season <- 0.5                           # if catch_formulation = discrete, 
                                            #       time at which fishing occurs
                                            #       0 = beginning, 1 = end
  }
  
  output = list(max_age, M, rec_age, af, bf, am, bm, a1f, L1f, a2f, L2f, Kf,
                a1m, L1m, a2m, L2m, Km, L50, k_mat, ldp, h, phi, sigma_R,
                rho_R, AMP, D, Fb, r, x, sp, c, b, fleets, alpha, beta, start, 
                F_fin, L50_up, L50_down, cf, switch, full, catch_form, season)
  
  return(output)
  
}