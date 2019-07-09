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
source("./code/R/stable_age_distribution.R")
source("./code/R/catch_at_age.R")
source("./code/R/effort_allocation.R")

##### Load life history characteristics for species ############################

par <- parameters("black rockfish")

max_age  <- par[[1]]                      # maximum age
M        <- par[[2]]                      # natural mortality
rec_age  <- par[[3]]                      # age at recruitment
af <- par[[4]];   bf  <- par[[5]]         # weight at length parameters (f)
am <- par[[6]];   bm  <- par[[7]]         # weight at length parameters (m)
a1f <- par[[8]];  L1f <- par[[9]];        # growth parameters (f)
a2f <- par[[10]]; L2f <- par[[11]]; 
Kf <- par[[12]]  
a1m <- par[[13]]; L1m <- par[[14]];       # growth parameters (m)
a2m <- par[[15]]; L2m <- par[[16]]; 
Km       <- par[[17]]  
L50      <- par[[18]]                     # length at 50% maturity
k_mat    <- par[[19]]                     # slope of maturity curve
ldp      <- par[[20]]                     # larval drift proportion
h        <- par[[21]]                     # steepness
phi      <- par[[22]]                     # unfished recruits per spawner
sigma_R  <- par[[23]]                     # recruitment standard deviation
rho_R    <- par[[24]]                     # recruitment autocorrelation
p        <- par[[25]]                     # adult movement proportion
D        <- par[[26]]                     # depletion
Fb       <- par[[27]]                     # fishing mortality to cause D
r        <- par[[28]]                     # proportion of positive transects 
                                          #       in PISCO monitoring data
x        <- par[[29]]                     # mean of positive transects
sp       <- par[[30]]                     # std of positive transects
c        <- par[[31]]                     # eggs produced per g, intercept
b        <- par[[32]]                     # eggs produced per g, slope

####### selectivity parameters #######
fleets   <- par[[33]]                     # fishery fleet names
alpha    <- par[[34]]                     # slope for upcurve
beta     <- par[[35]]                     # slope for downcurve
start    <- par[[36]]                     # length at initial vulnerability
F_fin    <- par[[37]]                     # F_fin for fishery, 0 if asymptotic
L50_up   <- par[[38]]                     # L50 for upcurve
L50_down <- par[[39]]                     # L50 for downcurve
cf       <- par[[40]]                     # fraction of fishery caught / fleet
switch   <- par[[41]]                     # length where selectivity switches 
                                          #       from upcurve to 1
full     <- par[[42]]                     # length at which downcurve starts


##### Population Dynamics - Non-Time Varying ###################################

# Set model parameters (fixed)
CR            <- 8                  # number of control rules
transects     <- 24                 # number of transects per PISCO protocol
                                    #     reserve implementation

# Set model parameters (flexible)
A             <- 5                  # number of areas, should be odd
time          <- 50                 # number of timesteps (years) before 
                                    #     reserve implementation
time2         <- 50                 # number of timesteps (years) after
allocation    <- 'equal'            # distribution of fishing effort (or 'IFD')
R0            <- 100000             # unfished recruitment, arbitrary value   
Init_size     <- 100000             # total population size at t = 1, 2

# Initialize arrays for time-varying dynamics
IA <- initialize_arrays(A, time, time2, R0, rec_age, max_age, L1f, L2f, Kf, a1f, 
                        a2f, af, bf, k_mat, Fb, L50, sigma_R, rho_R, fleets, 
                        alpha, beta, start, F_fin, L_50_up, L50_down, cf, 
                        switch, full, x, sp, M, CR, phi)

timeT            <- IA[[1]]       # total amount of timesteps (years)
E                <- IA[[2]]       # nominal fishing effort in each area 
age              <- IA[[3]]       # ages for which fish have recruited
n                <- IA[[4]]       # number of age classes
L                <- IA[[5]]       # Length at age, dim = 1*age
W                <- IA[[6]]       # Weight at age, dim = 1*age
Mat              <- IA[[7]]       # Fraction mature at age, dim = 1*age
m                <- IA[[8]]       # Age at which fraction mature > 0.5
S                <- IA[[9]]       # Selectivity at age
FM               <- IA[[10]]      # Fishing mortality rate, dim = age*area*time
N                <- IA[[11]]      # Population size, dim = age*area*time
SSB              <- IA[[12]]      # Spawning stock biomass, dim = area*time
abundance_all    <- IA[[13]]      # Abundance, dim = area*time
abundance_mature <- IA[[14]]      # Abundance, dim = area*time
biomass          <- IA[[15]]      # Biomass, dim = area*time
count_sp         <- IA[[16]]      # Species count when sampling, dim = area*time
nuS              <- IA[[17]]      # Sampling normal variable, dim = area*time*CR
Eps              <- IA[[18]]      # Epsilon vector, dim = area*time*CR
L0               <- IA[[19]]      # Length at age for stable age distribution
W0               <- IA[[20]]      # Weight at age for stable age distribution
catch            <- IA[[21]]      # Catch at age
yield            <- IA[[22]]      # Yield per area 
B0               <- IA[[23]]      # Unfished spawning stock biomass

##### Population Dynamics - Time Varying #######################################

for (t in 3:time) {
  
  for (a in 1:A) {
    
    # effort allocation
    E <- effort_allocation(a, t, allocation, A, E, biomass)
    
    # biology
    PD <- pop_dynamics(a, t, cr = 1, rec_age, max_age, n, SSB, N, W, Mat, A, R0, 
                       h, B0, Eps, sigma_R, Fb, E, S, M)
    SSB <- PD[[1]]
    FM  <- PD[[2]]
    N   <- PD[[3]]
    
    abundance_all[a, t, 1] <- sum(N[, a, t, 1])
    abundance_mature[a, t, 1] <- sum(N[m:(max_age-1), a, t, 1])
    biomass[a, t, 1] <- sum(N[, a, t, 1] * W)
    
    # sampling
    if (t > (time - 3)) {
      count_sp <- sampling(a, t + 3 - time, cr = 1, r, D, abundance_all, 
                           abundance_mature, transects, x, count_sp, nuS)
    }
    
    # fishing
    catch[, a, t, 1] <- catch_at_age(a, t, cr = 1, N, FM)
    N[, a, t, 1] <- N[, a, t, 1] - catch[, a, t, 1]
    yield[a, t, 1] <- sum(catch[, a, t, 1]*W)
    
  }
  
}

##### Implement Reserve, and apply control rules ###############################

for (cr in 1:CR) {
  
  for (a in 1:A) {
    
    for (t in 1:time2) {
      # effort allocation
      E <- effort_allocation(a, t, allocation, A, E, biomass)
      
      # biology
      PD <- pop_dynamics(a, t + time, cr, rec_age, max_age, n, SSB, N, W, Mat, A, 
                         R0, h, B0, Eps, sigma_R, Fb, E, S, M)
      SSB <- PD[[1]]
      FM  <- PD[[2]]
      N   <- PD[[3]]
      
      abundance_all[a, t + time, cr] <- sum(N[, a, t + time, cr])
      abundance_mature[a, t + time, cr] <- sum(N[m:(max_age-1), a, t + time, cr])
      
      biomass[a, t + time, cr] <- sum(N[, a, t + time, cr] * W)
      
      # sampling
      count_sp <- sampling(a, t + 3, cr, r, D, abundance_all, 
                           abundance_mature, transects, x, count_sp, nuS)
      
      # management
      E <- control_rule(a, t + 3, E, count_sp, x)
      
      # fishing
      catch[, a, t + time, r] <- catch_at_age(a, t + time, cr, N, FM)
      N[, a, t + time, cr] <- N[, a, t + time, cr] - catch[, a, t + time, cr]
      yield[a, t + time, cr] <- sum(catch[, a, t + time, cr]*W)
      
    }
    
  }
  
#### plot abundance, biomass, and yield over time for each area, once per CR ###
  
  a <- 1
  par(mfrow = c(1, 2))
  
  main_title <- sprintf("Control Rule %i", cr)
  
  # plot abundance (1000s of individuals) in blue
  plot(1:timeT, abundance_all[a, , cr]/1000, pch = 16, col = "deepskyblue3", 
       xlab = 'Time (years)', ylab = 'Abundance (1000s of individuals)',
       yaxt = 'n', ylim = c(0, 1500), xaxt = 'n', main = main_title)
  axis(1, seq(0, 100, 50))
  axis(2, seq(0, 1500, 500))
  
  # add red line for biomass (metric tons)
  lines(1:timeT, biomass[a, , cr]/1000, type = 'l', lwd = 2, col = "firebrick3")
  box()
  
  # plot yield over time (metric tons)
  plot(1:timeT, yield[a, , cr]/1000, type = 'l', lwd = 2, col = "forestgreen",
       xlab = 'Time (years)', ylab = 'Yield (metric tons)', 
       yaxt = 'n', ylim = c(0, 20), xaxt = 'n', main = main_title)
  axis(1, seq(0, 100, 50))
  axis(2, seq(0, 20, 5))
  box() 
  
}
