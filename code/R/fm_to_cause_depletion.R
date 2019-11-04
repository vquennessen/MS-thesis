# Calculate Fm to cause observed depletion

# clear environment
rm(list = ls())

species <- 'black rockfish 2003'
eq_time <- 150
true_dep <- 0.488
R0 <- 1e+5
stochasticity <- F

# source required functions
source("./parameters.R")
source("./length_at_age.R")
source("./weight_at_age.R")
source("./fraction_mature_at_age.R")
source("./fishing_mortality.R")
source("./epsilon.R")
source("./selectivity_at_age.R")
source("./spawning_stock_biomass.R")
source("./recruitment.R")
source("./pop_dynamics.R")
source("./initialize_arrays.R")
source("./sampling.R")
source("./density_ratio.R")
source("./management.R")
source("./control_rule.R")
source("./Leslie_SAD.R")
source("./catch_at_age.R")
source("./effort_allocation.R")
source("./equilibrium_SAD.R")
source("./movement.R")

# load species parameters
par <- parameters(species)

max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
rec_age                <- par[[3]]        # age at recruitment
af  <- par[[4]];   bf  <- par[[5]]        # weight at length parameters (f)
a1f <- par[[6]];  L1f  <- par[[7]]        # growth parameters (f)
a2f <- par[[8]];  L2f  <- par[[9]] 
Kf  <- par[[10]]  
L50                    <- par[[11]]       # length at 50% maturity
k_mat                  <- par[[12]]       # slope of maturity curve
ldp                    <- par[[13]]       # larval drift proportion
h                      <- par[[14]]       # steepness
phi                    <- par[[15]]       # unfished recruits per spawner
sigma_R                <- par[[16]]       # recruitment standard deviation
rho_R                  <- par[[17]]       # recruitment autocorrelation
AMP                    <- par[[18]]       # adult movement proportion
D                      <- par[[19]]       # depletion
Fb                     <- par[[20]]       # fishing mortality to cause D
r                      <- par[[21]]       # proportion of positive transects 
#       in PISCO monitoring data
x                      <- par[[22]]       # mean of positive transects
sp                     <- par[[23]]       # std of positive transects
c                      <- par[[24]]       # eggs produced per g, intercept
b                      <- par[[25]]       # eggs produced per g, slope

####### selectivity parameters #######
fleets                 <- par[[26]]       # fishery fleet names
alpha                  <- par[[27]]       # slope for upcurve
beta                   <- par[[28]]       # slope for downcurve
F_fin                  <- par[[29]]       # F_fin for fishery, 0 if asymptotic
L50_up                 <- par[[30]]       # L50 for upcurve
L50_down               <- par[[31]]       # L50 for downcurve
cf                     <- par[[32]]       # fraction of fishery caught / fleet

# Calculated values
age <- rec_age:max_age                          # applicable ages
n <- length(age)                                # number of age bins
L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f) # length at age
W <- weight_at_age(L, af, bf)                   # weight at age
Mat <- fraction_mature_at_age(n, k_mat, L, L50) # maturity at age
m <- age[min(which(Mat > 0.5))]                 # age at 50% mature
B0 <- R0/phi                                    # unfished spawning stock biomass
S <- selectivity_at_age(fleets, L, max_age, rec_age, alpha, L50_up, 
                        L50_down, F_fin, beta, n, cf, age)
Fb <- 0

# Initialize population size and catch arrays
# Dimensions = age * 1 * time * 1
N2 <- catch2 <- array(rep(0, n*eq_time), c(n, eq_time, 1, 1))

# Initialize biomass, SSB, and recruitment error
# Dimensions = 1 * time * 1
SSB2 <- Eps2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))

# Recruitment normal variable
# Dimensions = area * timeT * CR
if (stochasticity == T) {
  nuR2 <- array(rnorm(eq_time, 0, sigma_R), c(1, eq_time, 1, 1))
} else if (stochasticity == F) {
  nuR2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
}

# eps[1]
Eps2[1, 1, 1, 1] <- nuR2[1, 1, 1, 1]*sqrt(1 + rho_R^2)

# fill in rest of epsilon vector
for (t in 2:eq_time) {
  Eps2[1, t, 1, 1] <- rho_R*Eps2[1, t-1, 1, 1] + 
    nuR2[1, t, 1, 1]*sqrt(1 + rho_R^2)
}

# Fishing effort stays constant
E2 <- array(rep(1, eq_time), c(1, eq_time, 1, 1))

# Initialize FM and depletion levels
FM_values <- seq(from = 0, to = 1, by = 0.01)
fn <- length(FM_values)
dep <- rep(0, fn)

biomass2 <- array(rep(NA, fn*eq_time), c(fn, eq_time))

# Enter FM, N, abundance, and biomasses for time = 1 to rec_age
for (t in 1:rec_age) {
  N2[, a, t, 1, 1] <- SAD
  biomass2[a, t, 1, 1] <- sum(N2[, a, t, 1, 1] * W)
  catch2[a, t, 1, 1] <- 0
  SSB2[t, 1, 1] <- sum(N2[, a, t - rec_age, 1, 1]*W*Mat)
}

# Substitute in values for Fb to get depletion level
for (i in 1:fn) { 
  
  FM2 <- array(rep(FM_values[i], n*eq_time), c(n, 1, eq_time, 1, 1))
  
  # Step population forward in time with set fishing level
  for (t in (rec_age + 1):eq_time) {
    
    PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, rec_age, max_age, n, SSB2, 
                       N2, W, Mat, A, R0, h, B0, Eps2, sigma_R, Fb, E2, S, 
                       NM, FM2, m, abundance_all2, abundance_mature2, 
                       biomass2, fishing = F, nat_mortality)
    
    N2[, a, t, cr, nm]                <- PD[[2]]
    biomass2[a, t, cr, nm]            <- PD[[5]]

  }
  
}

dep <- 1 - (biomass2[eq_time, 1, 1, 1] / biomass2[eq_time, 1, 1, 1])

closest_FM <- FM_values[which.min(abs(dep - true_dep))] 

plot(FM_values, dep)
abline(v = closest_FM, col = 'red')
abline(h = true_dep, col = 'green')
print(closest_FM)
