# Calculate Fm to cause observed depletion

# clear environment
rm(list = ls())

species <- 'black rockfish 2003'
eq_time <- 150
true_dep <- 0.488
R0 <- 1e+5
stochasticity <- T

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
Fb <- 0

# Initialize population size and catch arrays
# Dimensions = age * 1 * time * 1
N2 <- catch2 <- array(rep(0, n*eq_time), c(n, eq_time))

# Initialize biomass, SSB, and recruitment error
# Dimensions = 1 * time * 1
SSB2 <- Eps2 <- array(rep(0, eq_time), c(eq_time))

# Recruitment normal variable
# Dimensions = area * timeT * CR
if (stochasticity == T) {
  nuR2 <- array(rnorm(eq_time, 0, sigma_R), c(eq_time))
} else if (stochasticity == F) {
  nuR2 <- array(rep(0, eq_time), c(eq_time))
}

# eps[1]
Eps2[1] <- nuR2[1]*sqrt(1 + rho_R^2)

# fill in rest of epsilon vector
for (t in 2:eq_time) {
  Eps2[t] <- rho_R*Eps2[t-1] + nuR2[t]*sqrt(1 + rho_R^2)
}

# Fishing effort stays constant
E2 <- array(rep(1, eq_time), c(1, eq_time))

# Initialize FM and depletion levels
FM_values <- seq(from = 0, to = 1, by = 0.001)
fn <- length(FM_values)
dep <- rep(0, fn)

biomass2 <- array(rep(NA, fn*eq_time), c(fn, eq_time))

# Enter FM, N, abundance, and biomasses for time = 1 to rec_age
for (t in 1:rec_age) {
  N2[, t] <- rep(10, n)
  biomass2[, t] <- sum(N2[, t] * W)
  catch2[, t] <- 0
  SSB2[t] <- sum(N2[, t - rec_age]*W*Mat)
}

# Substitute in values for Fb to get depletion level
for (i in 1:fn) { 
  
  # Step population forward in time with set fishing level
  for (t in (rec_age + 1):eq_time) {
    
    # set constant fishing mortality
    FM2 <- FM_values[i]
    
    # Calculate spawning stock biomass
    SSB2[t] <- sum(N2[, t - rec_age]*W*Mat)
    
    ##### Step population foward in time
    
    # Calculate recruitment and add recruits to population
    R1 <- (0.8 * R0 * h * SSB2[t-1]) / (0.2 * B0 * (1 - h) + (h - 0.2) * SSB2[t-1]) 
    N2[1, t] <- R1 * (exp(Eps2[t] - sigma_R^2 / 2))
    
    # Ages rec_age + 1 to max_age - 1
    for (j in 2:(n - 1)) {
      N2[j, t] <- N2[j - 1, t - 1] * exp(-1 * (FM2 + M))
    }
    
    # Final age bin
    N2[n, t] <- N2[n - 1, t - 1] * exp(-1 * (FM2 + M)) + 
      N2[n, t - 1] * exp(-1 * (FM2 + M)) 
    
    # fishing
    coeff <- FM2/(M + FM2)
    catch2[ , t] <- coeff * N2[ , t] * exp(-1*(M + FM2))
    
    N2[, t] <- N2[, t] - catch2[, t]
    biomass2[i, t] <- sum(N2[, t] * W)
    
  }
  
 # biomass3[i]<-biomass2[]
  # plot(1:eq_time, biomass2, main = i)
  
}

dep <- 1 - (biomass2[, eq_time] / biomass2[1, eq_time])

closest_FM <- FM_values[which.min(abs(dep - true_dep))] 

plot(FM_values, dep)
abline(v = closest_FM, col = 'red')
abline(h = true_dep, col = 'green')
print(closest_FM)