# Calculate Fm to cause observed depletion

# clear environment
rm(list = ls())

species <- 'black rockfish 2003'
eq_time <- 151
true_dep <- 0.488
R0 <- 1e+5
stochasticity <- F

# source required functions
source("./code/R/parameters.R")
source("./code/R/length_at_age.R")
source("./code/R/weight_at_age.R")
source("./code/R/fraction_mature_at_age.R")
source("./code/R/old_selectivity_at_age.R")
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
source("./code/R/Leslie_SAD.R")
source("./code/R/catch_at_age.R")
source("./code/R/effort_allocation.R")
source("./code/R/initial_size.R")
source("./code/R/vulnerability_to_gear.R")
source("./code/R/equilibrium_SAD.R")
source("./code/R/movement.R")

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
ldp                    <- par[[20]]       # larval drift proportion
h                      <- par[[21]]       # steepness
phi                    <- par[[22]]       # unfished recruits per spawner
sigma_R                <- par[[23]]       # recruitment standard deviation
rho_R                  <- par[[24]]       # recruitment autocorrelation
AMP                    <- par[[25]]       # adult movement proportion
D                      <- par[[26]]       # depletion
r                      <- par[[28]]       # proportion of positive transects 
                                          #       in PISCO monitoring data
x                      <- par[[29]]       # mean of positive transects
sp                     <- par[[30]]       # std of positive transects
c                      <- par[[31]]       # eggs produced per g, intercept
b                      <- par[[32]]       # eggs produced per g, slope

####### selectivity parameters #######
fleets                 <- par[[33]]       # fishery fleet names
alpha                  <- par[[34]]       # slope for upcurve
beta                   <- par[[35]]       # slope for downcurve
start                  <- par[[36]]       # length at initial vulnerability
F_fin                  <- par[[37]]       # F_fin for fishery, 0 if asymptotic
L50_up                 <- par[[38]]       # L50 for upcurve
L50_down               <- par[[39]]       # L50 for downcurve
cf                     <- par[[40]]       # fraction of fishery caught / fleet
switch                 <- par[[41]]       # length where selectivity switches 
                                          #       from upcurve to 1
full                   <- par[[42]]       # length at which downcurve starts
catch_form             <- par[[43]]       # discrete or continuous catch
season                 <- par[[44]]       # if catch_formulation = discrete, 
                                          #       time at which fishing occurs:
                                          #       0 at start, 1 at end of year

# Calculated values
age <- rec_age:max_age                          # applicable ages
n <- length(age)                                # number of age bins
L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f) # length at age
W <- weight_at_age(L, af, bf)                   # weight at age
Mat <- fraction_mature_at_age(n, k_mat, L, L50) # maturity at age
m <- age[min(which(Mat > 0.5))]                 # age at 50% mature
B0 <- R0/phi                                    # unfished spawning stock biomass
S <- old_selectivity_at_age(L, fleets, alpha, beta, # selectivity at age
                            start, F_fin, L50_up, L50_down, cf, switch, full)

# Initialize population size and catch arrays
# Dimensions = age * 1 * time * 1
N2 <- catch2 <- FM2 <- array(rep(0, n*eq_time), c(n, 1, eq_time, 1))

# Initialize biomass, SSB, and recruitment normal variable arrays
# Dimensions = 1 * time * 1
biomass2 <- SSB2 <- E2 <- array(rep(0, eq_time), c(1, eq_time, 1))
abundance_all2 <- abundance_mature2 <- array(rep(0, eq_time), c(1, eq_time, 1))

# Recruitment normal variable
# Dimensions = area * timeT * CR
if (stochasticity == T) {
  nuR2 <- array(rnorm(1*eq_time*1, 0, sigma_R), c(1, eq_time, 1))
} else if (stochasticity == F) {
  nuR2 <- array(rep(0, 1*eq_time*1), c(1, eq_time, 1))
}

# Recruitment error
# Dimensions = area * eq_time - 1 * CR
Eps2 <- epsilon(A = 1, eq_time - 1, CR = 1, nuR2, rho_R)

# Fishing effort stays constant
E2[, 1:eq_time, ] <- rep(1, eq_time)

# Initialize FM and depletion levels
fb_values <- seq(from = 0, to = 1, by = 0.01)
fn <- length(fb_values)
dep <- rep(0, fn)

# Enter FM, N, abundance, and biomasses for time = 1 to rec_age
for (t in 1:rec_age) {
  FM2[, 1, t, 1] <- 0
  N2[, 1, t, 1] <- rep(1, n)
  biomass2[1, t, 1] <- sum(N2[, 1, t, 1] * W)
  catch2[, 1, t, 1] <- 0
  SSB2[1, t, 1] <- spawning_stock_biomass(a = 1, t, cr = 1, rec_age, N2, 
                                          W, Mat)
}

# Substitute in values for Fb to get depletion level
for (i in 1:fn) { 
  
  # # Start each age class with Stable Age Distribution
  # # Enter FM, N, abundance, and biomasses for time = 1 to rec_age + 1
  # # Dimensions = age * area * time * CR
  # SAD <- equilibrium_SAD(1, 1, A = 1, rec_age, max_age, n, W, R0,
  #                        Mat, h, B0, Eps2, sigma_R, Fb = fb_values[i], 
  #                        S, M, season, catch_form, eq_time, m, stochasticity, 
  #                        rho_R)
  # 
  # # Enter FM, N, abundance, and biomasses for time = 1 to rec_age
  #   for (t in 1:rec_age) {
  #       FM2[, 1, t, 1] <- fishing_mortality(a = 1, t, cr = 1, FM2, A = 1, 
  #                                           Fb = fb_values[i], E2, S)
  #       N2[, 1, t, 1] <- SAD
  #       biomass2[1, t, 1] <- sum(N2[, 1, t, 1] * W)
  #       catch2[, 1, t, 1] <- catch_at_age(a = 1, t, cr = 1, FM2, M, N2, A = 1, 
  #                                         Fb = fb_values[i], E2, catch2, 
  #                                         catch_form, season)
  #       N2[, 1, t, 1] <- N2[, 1, t, 1] - catch2[, 1, t, 1]
  #       SSB2[1, t, 1] <- spawning_stock_biomass(a = 1, t, cr = 1, rec_age, N2, 
  #                                               W, Mat)
  #     }

  # Step population forward in time with set fishing level
  for (t in (rec_age + 1):(eq_time - 1)) {
    
    # biology
    PD <- pop_dynamics(a = 1, t, cr = 1, rec_age, max_age, n, SSB2, N2, W, 
                       Mat, A = 1, R0, h, B0, Eps2, sigma_R, 
                       Fb = fb_values[i], E2, S, M, FM2, m, abundance_all2, 
                       abundance_mature2, biomass2)
    
    SSB2                <- PD[[1]]
    FM2                 <- PD[[2]]
    N2                  <- PD[[3]]
    biomass2            <- PD[[6]]
    
    # fishing
    catch2[, 1, t, 1] <- catch_at_age(a = 1, t, cr = 1, FM2, M, N2, A = 1, 
                                      Fb = fb_values[i], E2, catch2, 
                                      catch_form, season)
    N2[, 1, t, 1] <- N2[, 1, t, 1] - catch2[, 1, t, 1]
    biomass2[1, t, 1] <- sum(N2[, 1, t, 1] * W)
    
  }
  
  dep[i] <- 1 - (biomass2[eq_time - 1] / B0)
  
}

closest_Fb <- fb_values[which.min(abs(dep - true_dep))] 

plot(fb_values, dep)
abline(v = closest_Fb, col = 'red')
print(closest_Fb)