selectivity_values <- function(species, fleet, gender, Sa, La, L50up, L50down) {
  
  # source required functions
  source("./code/R/parameters.R")
  source("./code/R/length_at_age.R")

  species <- 'black rockfish 2015'
  fleet <- 'trawl'
  gender <- 'female'
  Sa <- 0.3
  La <- 6
  L50up <- 7
  
  # load species parameters
  par <- parameters(species)
  
  max_age                <- par[[1]]        # maximum age
  M                      <- par[[2]]        # natural mortality
  rec_age                <- par[[3]]        # age at recruitment
  af  <- par[[4]];   bf  <- par[[5]]        # weight at length parameters (f)
  am  <- par[[6]];   bm  <- par[[7]]        # weight at length parameters (m)
  a1f <- par[[8]];  L1f  <- par[[9]]        # growth parameters (f)
  a2f <- par[[10]]; L2f  <- par[[11]] 
  Kf  <- par[[12]]  
  a1m <- par[[13]]; L1m  <- par[[14]]       # growth parameters (m)
  a2m <- par[[15]]; L2m  <- par[[16]]  
  Km                     <- par[[17]]  
  L50                    <- par[[18]]       # length at 50% maturity
  k_mat                  <- par[[19]]       # slope of maturity curve
 
  # Calculated values
  age <- rec_age:max_age                          # applicable ages
  n <- length(age)                                # number of age bins
  L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f) # length at age

  alpha <- (-1 * log(1/Sa - 1)) / (L[L50up] - La)
  
}