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
  
  if (species == "black rockfish") {
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
    k_mat <- 0.4103                         # slope of maturity curve
    ldp <- 0.1                              # larval drift proportion
    R0 <- 100000                            # unfished recruitment
    h <- 0.65                               # steepness
    phi <- 1.1                              # unfished recruits per spawner
    sigma_R <- 0.5                          # recruitment standard deviation
    rho_R <- 0                              # recruitment autocorrelation
    p <- 0.1                                # adult movement proportion
    D <- 0.488                              # depletion
    Fb <- 0.2                               # fishing mortality to cause D
    r <- 0.77                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    x <- 15.42                              # mean of positive transects
    sp <- 16.97                             # std of positive transects
                                            
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
                                            
  }
  
  output = list(max_age, M, rec_age, af, bf, am, bm, a1f, L1f, a2f, L2f, Kf,
                a1m, L1m, a2m, L2m, Km, L50, k_mat, ldp, R0, h, phi, sigma_R,
                rho_R, p, D, Fb, r, x, sp, fleets, alpha, beta, start, F_fin, 
                L50_up, L50_down, cf, switch, full)
  
  return(output)
  
}