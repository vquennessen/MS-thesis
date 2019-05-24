#' Runs age-structured population dynamics biological sub-model. 

##### Source functions #########################################################

source("./code/R/parameters.R")
source("./code/R/length_at_age.R")
source("./code/R/weight_at_age.R")
source("./code/R/fraction_mature_at_age.R")
source("./code/R/selectivity_at_age.R")
source("./code/R/fishing_mortality.R")
source("./code/R/epsilon.R")
source("./code/R/spawning_stock_biomass.R")
source("./code/R/recruitment.R")
source("./code/R/pop_dynamics.R")

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
B0 <- par[[32]]                           # unfished spawning stock biomass, in
                                          #       metric tons

                                          ####### selectivity parameters #######
fleets <- par[[33]]                       # fishery fleet names
alpha <- par[[34]]                        # slope for upcurve
beta <- par[[35]]                         # slope for downcurve
start <- par[[36]]                        # length at initial vulnerability
F_fin <- par[[37]]                        # F_fin for fishery, 0 if asymptotic
L50_up <- par[[38]]                       # L50 for upcurve
L50_down <- par[[39]]                     # L50 for downcurve
cf <- par[[40]]                           # fraction of fishery caught / fleet
switch <- par[[41]]                       # length where selectivity switches 
                                          #       from upcurve to 1
full <- par[[42]]                         # length at which downcurve starts


##### Population Dynamics - Non-Time Varying ###################################

A <- 5                               # number of areas
time <- 50                           # number of timesteps (years)
E <- 0.10                            # nominal fishing effort in each area 
age <- rec_age:max_age               # ages for which fish have recruited
n <- length(age)

# Length at age
# Dimensions = 1 * age
L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f)

# Weight at age
# Dimensions = 1 * age
W <- weight_at_age(L, af, bf)

# Maturity at age
# Dimensions = 1 * age
M <- fraction_mature_at_age(n, k_mat, L, L50)

# Selectivity at age
# Dimensions = 1 * age
S <- selectivity_at_age(L, fleets, alpha, beta, start, F_fin, L50_up, L50_down, 
                        cf, switch, full)

# Fishing mortality
# Initialize array
# Dimensions = age * area * time
FM <- array(rep(0, n*A*time), c(n, A, time))

# First year of fishing mortality
FM[, , 1] <- fishing_mortality(A, Fb, E, S)

# Recruitment error
# Dimensions = 1 * time
e <- epsilon(time, sigma_R, rho_R)

# Initialize age-structured population size matrix
# Dimensions = age * area * time
N <- array(rep(0, n*A*time), c(n, A, time))

# Initial age structure
s <- 100          # start with 100 individuals in each area at t = 1
N[, , 1] <- array(rep(s, n), c(1, 1, n))

# Initialize spawning stock biomass array
# Dimensions = area * time
B <- array(rep(NA, A*time), c(A, time))

# Initialize recruitment vector, 
# Dimensions = area * time
R <- array(rep(0, A*time), c(A, time))

##### Population Dynamics - Time Varying #######################################

for (a in 1:A) {
  
  for (t in 2:time) {
    
    PD <- pop_dynamics(a, t, B, N, W, M, A, R0, h, B0, e, sigma_R, Fb, E, S)
    
  }
  
}
