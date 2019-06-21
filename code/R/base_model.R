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
c <- par[[33]]                            # eggs produced per g, intercept
b <- par[[34]]                            # eggs produced per g, slope

####### selectivity parameters #######
fleets <- par[[35]]                       # fishery fleet names
alpha <- par[[36]]                        # slope for upcurve
beta <- par[[37]]                         # slope for downcurve
start <- par[[38]]                        # length at initial vulnerability
F_fin <- par[[39]]                        # F_fin for fishery, 0 if asymptotic
L50_up <- par[[40]]                       # L50 for upcurve
L50_down <- par[[41]]                     # L50 for downcurve
cf <- par[[42]]                           # fraction of fishery caught / fleet
switch <- par[[43]]                       # length where selectivity switches 
                                          #       from upcurve to 1
full <- par[[44]]                         # length at which downcurve starts


##### Population Dynamics - Non-Time Varying ###################################

# Set model parameters
A             <- 5                  # number of areas
time          <- 50                 # number of timesteps (years) before 
                                    #     reserve implementation
time2         <- 50                 # number of timesteps (years) after
transects     <- 24                 # number of transects per PISCO protocol
                                    #     reserve implementation
init_effort   <- 1               # nominal fishing effort in each area
initial       <- 1000000            # total population size at t = 1, 2
CR            <- 8                  # number of control rules
allocation    <- 'equal'            # distribution of fishing effort

# Initialize arrays for time-varying dynamics
IA <- initialize_arrays(time, time2, init_effort, rec_age, max_age, L1f, 
                        L2f, Kf, a1f, a2f, af, bf, k_mat, Fb, L50, 
                        sigma_R, rho_R, fleets, alpha, beta, start, F_fin, 
                        L_50_up, L50_down, cf, switch, full, A, x, sp, 
                        initial, M, CR)

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
e                <- IA[[11]]      # Recruitment error, dim = 1*time
N                <- IA[[12]]      # Population size, dim = age*area*time
SSB              <- IA[[13]]      # Spawning stock biomass, dim = area*time
R                <- IA[[14]]      # Recruitment, dim = area*time
abundance_all    <- IA[[15]]      # Abundance, dim = area*time
abundance_mature <- IA[[16]]      # Abundance, dim = area*time
biomass          <- IA[[17]]      # Biomass, dim = area*time
count_sp         <- IA[[18]]      # Species count when sampling, dim = area*time
nu               <- IA[[19]]      # Sampling normal variable, dim = area*time
L0               <- IA[[20]]      # Length at age for stable age distribution
W0               <- IA[[21]]      # Weight at age for stable age distribution
catch            <- IA[[22]]      # Catch at age
yield            <- IA[[23]]      # Yield per area 

##### Population Dynamics - Time Varying #######################################

for (t in 3:time) {
  
  for (a in 1:A) {
    
    # effort allocation
    E <- effort_allocation(allocation, E, biomass, a, t)
    
    # biology
    PD <- pop_dynamics(a, t, y = 1, rec_age, max_age, n, SSB, N, W, Mat, A, R0, 
                       h, B0, e, sigma_R, Fb, E, S, M)
    SSB <- PD[[1]]
    R   <- PD[[2]]
    FM  <- PD[[3]]
    N   <- PD[[4]]
    
    abundance_all[a, t] <- sum(N[, a, t])
    abundance_mature[a, t] <- sum(N[m:(max_age-1), a, t])
    biomass[a, t] <- sum(N[, a, t] * W)
    
    # sampling
    if (t > (time - 3)) {
      count_sp <- sampling(a, t + 3 - time, y = 1, r, D, abundance_all, 
                           abundance_mature, transects, x, count_sp, nu)
    }
    
    # fishing
    catch[, a, t] <- catch_at_age(a, t, N, W, FM)
    N[, a, t] <- N[, a, t] - catch[, a, t]
    yield[a, t] <- sum(catch[, a, t])
    
  }
  
}

##### Implement Reserve, and apply control rules ###############################

for (y in 1:CR) {
  
  for (a in 1:A) {
    
    for (t in 1:time2) {
      # effort allocation
      E <- effort_allocation(allocation, E, biomass, a, t)
      
      # biology
      PD <- pop_dynamics(a, t + time, y, rec_age, max_age, n, SSB, N, W, Mat, A, 
                         R0, h, B0, e, sigma_R, Fb, E, S, M)
      SSB <- PD[[1]]
      R   <- PD[[2]]
      FM  <- PD[[3]]
      N   <- PD[[4]]
      
      abundance_all[a, t + time] <- sum(N[, a, t + time])
      abundance_mature[a, t + time] <- sum(N[m:(max_age-1), a, t + time])
      
      biomass[a, t + time] <- sum(N[, a, t + time] * W)
      
      # sampling
      count_sp <- sampling(a, t + 3, y, r, D, abundance_all, 
                           abundance_mature, transects, x, count_sp, nu)
      
      # management
      E <- control_rule(a, t + 3, E, count_sp, x)
      
      # fishing
      catch[, a, t + time] <- catch_at_age(a, t + time, N, W, FM)
      N[, a, t + time] <- N[, a, t + time] - catch[, a, t + time]
      yield[a, t + time] <- sum(catch[, a, t + time])
      
    }
    
  }
  
  # plot abundance, biomass, and yield over time for each area, once per CR
  
  # par(mfrow = c(5, 2))
  
  for (a in 1:A) {
    
    #par(mfrow = c(2, 1))
    
    par(mfrow = c(1, 2))
       
    main_title <- sprintf("Control Rule %i, Area %i", y, a)
    
    # plot abundance (1000s of individuals) in blue
    plot(1:timeT, abundance_all[a, ]/1000, pch = 16, col = "deepskyblue3", 
         xlab = 'Time (years)', ylab = 'Abundance (1000s of individuals)',
         yaxt = 'n', ylim = c(0, 1500), xaxt = 'n', main = main_title)
    axis(1, seq(0, 100, 50))
    axis(2, seq(0, 1500, 500))
    
    # add red line for biomass (metric tons)
    lines(1:timeT, biomass[a, ]/1000, type = 'l', lwd = 2, col = "firebrick3")
    box()
    
    # plot yield over time (metric tons)
    plot(1:timeT, yield[a, ]/1000, type = 'l', lwd = 2, col = "forestgreen",
         xlab = 'Time (years)', ylab = 'Yield (metric tons)', 
         yaxt = 'n', ylim = c(0, 20), xaxt = 'n')
    axis(1, seq(0, 100, 50))
    axis(2, seq(0, 20, 5))
    box() 
    
  }
  
}