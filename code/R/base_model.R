#' Runs age-structured population dynamics biological sub-model. 

##### Source functions #########################################################

source("./code/R/parameters.R")
source("./code/R/length_at_age.R")
source("./code/R/weight_at_age.R")
source("./code/R/fraction_mature_at_age.R")

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
alpha <- par[[32]]                        # slope for upcurve
beta <- par[[33]]                         # slope for downcurve
cf <- par[[34]]                           # fraction of fishery caught / fleet
F_fin <- par[[35]]                        # F_fin for fishery, 0 if asymptotic
L50_up <- par[[36]]                       # L50 for upcurve
L50_down <- par[[37]]                     # L50 for downcurve

##### Population Dynamics ######################################################

s <- 5            # number of areas
t <- 50           # number of timesteps (years)
E <- 0.10         # nominal fishing effort in each area 
# TODO: update E value

# Length at age
L <- length_at_age(max_age, L1f, L2f, Kf, a1f, a2f)

# Weight at age
W <- weight_at_age(L, af, bf)

# Maturity at length
M <- fraction_mature_at_age(max_age, k_mat, L, L50)


# Based on Babcock & MacCall (2011): Eq. (4)
epsilon <- array(rep(0, t), t)
nu <- rnorm(t, 0, sigma_R)
for (i in 2:t) {
  epsilon[i] <- rho_R*epsilon[i-1] + nu[i]*sqrt(1 + rho_R^2)
}

# Spawning Stock Biomass
B <- array(rep(0, s*t), c(s, t))

# Recruitment
# Based on Babcock & MacCall (2011): Eq. (3)
R <- array(rep(0, s*t), c(s, t))

# Catchability
# Based on Babcock & MacCall (2011): Eq. (6)
q <- (s*Fb)/(s*E)

# Selectivity

# Fishing Mortality
# Based on Babcock & MacCall (2011): Eq. (5)
f <- array(rep(0, length(n)*s*t), c(length(n), s, t))

# Population size for different ages
# Based on Babcock & MacCall (2011): Eq. (1)
n <- rec_age:max_age
N <- array(rep(0, length(n)*s*t), c(length(n), s, t))
