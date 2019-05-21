#' Runs age-structured population dynamics biological sub-model. 

##### Source functions #########################################################

source("./code/R/parameters.R")
source("./code/R/length_at_age.R")
source("./code/R/weight_at_age.R")
source("./code/R/fraction_mature_at_age.R")
source("./code/R/selectivity_at_age.R")
source("./code/R/fishing_mortality.R")
source("./code/R/epsilon.R")
source("./code/R/recruitment.R")

##### Load life history characteristics for species ############################

par <- parameters("black rockfish")

max_age <- par[[1]]                       # maximum age
M <- par[[2]]                             # natural mortality
rec_age <- par[[3]]                       # age at recruitment
af <- par[[4]]; bf <- par[[5]]            # weight at length parameters (f)
am <- par[[6]]; bm <- par[[7]]            # weight at length parameters (m)
a1f <- par[[8]]; L1f <- par[[9]];         # growth parameters (f)
a2f <- par[[10]]; L2f <- par[[11]]; 
Kf <- par[[12]]  
a1m <- par[[13]]; L1m <- par[[14]];       # growth parameters (m)
a2m <- par[[15]]; L2m <- par[[16]]; 
Km <- par[[17]]  
L50 <- par[[18]]                          # length at 50% maturity
k_mat <- par[[19]]                        # slope of maturity curve
ldp <- par[[20]]                          # larval drift proportion
R0 <- par[[21]]                           # unfished recruitment
h <- par[[22]]                            # steepness
phi <- par[[23]]                          # unfished recruits per spawner
sigma_R <- par[[24]]                      # recruitment standard deviation
rho_R <- par[[25]]                        # recruitment autocorrelation
p <- par[[26]]                            # adult movement proportion
D <- par[[27]]                            # depletion
Fb <- par[[28]]                           # fishing mortality to cause D
r <- par[[29]]                            # Proportion of positive transects 
                                          #       in PISCO monitoring data
x <- par[[30]]                            # mean of positive transects
sp <- par[[31]]                           # std of positive transects

                                          ##### selectivity parameters #########
fleets <- par[[32]]                       # fishery fleet names
alpha <- par[[33]]                        # slope for upcurve
beta <- par[[34]]                         # slope for downcurve
start <- par[[35]]                        # length at initial vulnerability
F_fin <- par[[36]]                        # F_fin for fishery, 0 if asymptotic
L50_up <- par[[37]]                       # L50 for upcurve
L50_down <- par[[38]]                     # L50 for downcurve
cf <- par[[39]]                           # fraction of fishery caught / fleet
switch <- par[[40]]                       # length where selectivity switches 
                                          #       from upcurve to 1
full <- par[[41]]                         # length at which downcurve starts


##### Population Dynamics ######################################################

A <- 5            # number of areas
t <- 50           # number of timesteps (years)
E <- 0.10         # nominal fishing effort in each area 

# Length at age
L <- length_at_age(max_age, L1f, L2f, Kf, a1f, a2f)

# Weight at age
W <- weight_at_age(L, af, bf)

# Maturity at age
M <- fraction_mature_at_age(max_age, k_mat, L, L50)

# Selectivity at age
S <- selectivity_at_age(L, fleets, alpha, beta, start, F_fin, L50_up, L50_down, 
                        cf, switch, full)

# Fishing mortality
FM <- fishing_mortality(A, Fb, E, S)

# Recruitment error
e <- epsilon(t, sigma_R, rho_R)

# Initialize age-structured population size matrix
# Each row is a timestep, each column is an age
N <- array(rep(0, t*max_age), c(t, max_age))

# Initial age structure
n <- 100          # start with 100 individuals in each age
N[1, ] <- array(rep(n, max_age), c(1, max_age))

# Recruitment


# Population size for different ages
# Based on Babcock & MacCall (2011): Eq. (1)
n <- rec_age:max_age
N <- array(rep(0, length(n)*s*t), c(length(n), s, t))
