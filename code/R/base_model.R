#' Runs age-structured population dynamics biological sub-model. 

rm(list = ls())

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
source("./code/R/initialize_arrays.R")
source("./code/R/sampling.R")
source("./code/R/density_ratio.R")
source("./code/R/management.R")
source("./code/R/control_rule.R")

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
r <- par[[29]]                            # proportion of positive transects 
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

# Set model parameters
A         <- 5                       # number of areas
time      <- 50                      # number of timesteps (years) before 
                                     #     reserve implementation
time2     <- 50                      # number of timesteps (years) after
                                     #     reserve implementation
timeT     <- time + time2            # total amount of timesteps (years)
E         <- rep(0.10, A)            # nominal fishing effort in each area 
age       <- rec_age:max_age         # ages for which fish have recruited
n         <- length(age)             # number of age classes
transects <- 24                      # number of transects per PISCO protocol
initial   <- 1                       # number in each age class at t = 1, 2
CR        <- 8                       # number of control rules

# Initialize arrays for time-varying dynamics
IA <- initialize_arrays(L1f, L2f, Kf, a1f, a2f, af, bf, k_mat, Fb,
                        L50, sigma_R, rho_R, fleets, alpha, beta, start, 
                        F_fin, L_50_up, L50_down, cf, switch, full, age, 
                        n, A, time, time2, E, x, sp, initial)

L                <- IA[[1]]       # Length at age, dim = 1*age
W                <- IA[[2]]       # Weight at age, dim = 1*age
Mat              <- IA[[3]]       # Fraction mature at age, dim = 1*age
m                <- IA[[4]]       # Age at which fraction mature > 0.5
S                <- IA[[5]]       # Selectivity at age
FM               <- IA[[6]]       # Fishing mortality rate, dim = age*area*time
e                <- IA[[7]]       # Recruitment error, dim = 1*time
N                <- IA[[8]]       # Population size, dim = age*area*time
SSB              <- IA[[9]]       # Spawning stock biomass, dim = area*time
R                <- IA[[10]]      # Recruitment, dim = area*time
abundance_all    <- IA[[11]]      # Abundance, dim = area*time
abundance_mature <- IA[[12]]      # Abundance, dim = area*time
biomass          <- IA[[13]]      # Biomass, dim = area*time
count_sp         <- IA[[14]]      # Species count when sampling, dim = area*time
nu               <- IA[[15]]      # Sampling normal variable, dim = area*time

##### Population Dynamics - Time Varying #######################################
 
for (t in 3:time) {
  
  for (a in 1:A) {
    
    PD <- pop_dynamics(a, t, rec_age, max_age, n, SSB, N, W, Mat, A, R0, h, 
                       B0, e, sigma_R, Fb, E, S, M)
    SSB <- PD[[1]]
    R   <- PD[[2]]
    FM  <- PD[[3]]
    N   <- PD[[4]]
    
    abundance_all[a, t] <- sum(N[, a, t])
    abundance_mature[a, t] <- sum(N[m:(max_age-1), a, t])
    biomass[a, t] <- sum(N[, a, t] * W)
    
    if (t > (time - 3)) {
      count_sp <- sampling(a, t + 3 - time, r, D, abundance_all, 
                           abundance_mature, transects, x, count_sp, nu)
    }
    
  }
  
}

##### Implement Reserve, and apply control rules ###############################

for (x in 1:CR) {
  
  for (t in 1:time2) {
    
    for (a in 1:A) {
      
      PD <- pop_dynamics(a, t + time, rec_age, max_age, n, SSB, N, W, Mat, A, 
                         R0, h, B0, e, sigma_R, Fb, E, S, M)
      SSB <- PD[[1]]
      R   <- PD[[2]]
      FM  <- PD[[3]]
      N   <- PD[[4]]
      
      abundance_all[a, t] <- sum(N[, a, t])
      abundance_mature[a, t] <- sum(N[m:(max_age-1), a, t])
      
      biomass[a, t] <- sum(N[, a, t] * W)
      
      count_sp <- sampling(a, t + 3, r, D, abundance_all, 
                           abundance_mature, transects, x, count_sp, nu)
      
      E <- control_rule(a, t + 3, E, count_sp, x)
    }
    
  }
  
}