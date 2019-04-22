#' Loads life history characteristic parameters for selected species
#' 
#' @param species a character string
#' @param max_age maximum age, a number (> 0)
#' @param M natural mortality, a number (0, 1)
#' @param rec_age age at recruitment, a number (> 0, < max_age)
#' @param af female weight at length parameter, a number (> 0)
#' @param bf female weight at length parameter, a number (~ 3)
#' @param am male weight at length parameter, a number (> 0)
#' @param bm male weight at length parameter, a number (~ 3)
#' @param a1f female age 1, a number (> 0)
#' @param L1f female length at age 1, a number (> 0)
#' @param a2f female age 2, a number (> a1f)
#' @param L2f female length at age 2, a number (> L1f)
#' @param Kf female growth rate, a number (0, 1)
#' @param a1m male age 1, a number (> 0)
#' @param L1m male length at age 1, a number (> 0)
#' @param a2m male age 2, a number (> a1m)
#' @param L2m male length at age 2, a number (> L1m)
#' @param Km male growth rate, a number (0, 1)
#' @param L50 length at 50% maturity, a number (> 0)
#' @param k_mat slope of maturity curve, a number
#' @param ldp larval drift proportion, a number (0, 1)
#' @param R0 unfished recruitment, a number
#' @param h steepness, a number (0, 1)
#' @param phi unfished recruits per spawner, a number
#' @param sigma_R recruitment standard deviation, a number
#' @param rho_R recruitment autocorrelation, a number
#' @param p adult movement proportion, a number (0, 1)
#' @param D depletion, a number (0, 1)
#' @param Fb fishing mortality to cause D, a number (0, 1)
#' @param r proportion of positive transects in monitoring data, a number (0, 1)
#' @param x mean of positive transects, a number
#' @param sp standard deviation of positive transects, a number
#' 
#' @return the parameters needed to run the popDyn function

parameters = function(study_species) {
  
  # Error handling
  
  if (is.character(study_species) == 0) {
    stop("Study species must be a character string.")
  }
  
  if (study_species == "black rockfish") {
    # source: Cope et al., 2015
    max_age <- 40                           ## maximum age
    M <- 0.17                               ## natural mortality
    rec_age <- 3                            ## age at recruitment
    af <- 2.6e-5; bf <- 2.88                ## weight at length parameters (f)
    am <- 2.58e-5; bm <- 2.89               ## weight at length parameters (m)
    a1f <- 1; L1f <- 20.32;                 ## growth parameters (f)
    a2f <- 40; L2f <- 49.67; Kf <- 0.21  
    a1m <- 1; L1m <- 17.47;                 ## growth parameters (m)
    a2m <- 40; L2m <- 43.27; Km <- 0.34  
    L50 <- 43.69                            ## length at 50% maturity
    k_mat <- -0.66                          ## slope of maturity curve
    ldp <- 0.1                              # larval drift proportion
    R0 <- 3666                              ## unfished recruitment
    h <- 0.77                               ## steepness
    phi <- 1.1                              # unfished recruits per spawner
    sigma_R <- 0.5                          ## recruitment standard deviation
    rho_R <- 0                              # recruitment autocorrelation
    p <- 0.1                                # adult movement proportion
    D <- 0.604                              ## depletion
    Fb <- 0.2                               # fishing mortality to cause D
    r <- 0.77                               # Proportion of positive transects 
                                            #       in PISCO monitoring data
    x <- 15.42                              # mean of positive transects
    sp <- 16.97                             # std of positive transects
  }
  
}