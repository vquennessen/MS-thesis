# clear environment
rm(list = ls())

##### specify parameters for this run #####
species <- 'BR2003'
eq_time <- 150
true_dep <- 0.488
R0 <- 1e+5
stochasticity <- F
recruitment_mode <- 'pool'
A <- 5

##### source required functions #####
source("./parameters.R")
source("./length_at_age.R")
source("./weight_at_age.R")
source("./fraction_mature_at_age.R")
source("./selectivity_at_age.R")
source("./pop_dynamics.R")
source("./equilibrium_SAD.R")
source("./epsilon.R")
source("./recruitment.R")
source("./fishing_mortality.R")

##### load species parameters #####
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
h                      <- par[[14]]       # steepness
phi                    <- par[[15]]       # unfished recruits per spawner
sigma_R                <- par[[16]]       # recruitment standard deviation
rho_R                  <- par[[17]]       # recruitment autocorrelation
                                          #       in PISCO monitoring data
sp                     <- par[[23]]       # std of positive transects

####### selectivity parameters #######
fleets                 <- par[[26]]       # fishery fleet names
alpha                  <- par[[27]]       # slope for upcurve
beta                   <- par[[28]]       # slope for downcurve
F_fin                  <- par[[29]]       # F_fin for fishery, 0 if asymptotic
A50_up                 <- par[[30]]       # L50 for upcurve
A50_down               <- par[[31]]       # L50 for downcurve
cf                     <- par[[32]]       # fraction of fishery caught / fleet

##### Calculate set values #####
age <- rec_age:max_age                          # applicable ages
n <- length(age)                                # number of age bins
L <- length_at_age(rec_age, max_age, L1f, L2f, 
                   Kf, a1f, a2f, all_ages = F)  # length at age
W <- weight_at_age(L, af, bf)                   # weight at age
Mat <- fraction_mature_at_age(n, k_mat, L, L50) # maturity at age
m <- age[min(which(Mat > 0.5))]                 # age at 50% mature
B0 <- R0/phi                                    # unfished spawning stock biomass
S <- selectivity_at_age(fleets, L, max_age, rec_age, alpha, A50_up, 
                        A50_down, F_fin, beta, n, cf, age)

# Recruitment normal variable
# Dimensions = area * timeT * CR
if (stochasticity == T) {
  nuR2 <- array(rnorm(eq_time, 0, sigma_R), c(1, eq_time, 1, 1))
} else if (stochasticity == F) {
  nuR2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
}

# initialize epsilon vector
Eps2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
# eps[1]
Eps2[1, 1, 1, 1] <- nuR2[1, 1, 1, 1]*sqrt(1 + rho_R^2)
# fill in rest of epsilon vector
for (t in 2:eq_time) {
  Eps2[1, t, 1, 1] <- rho_R*Eps2[1, t-1, 1, 1] + 
    nuR2[1, t, 1, 1]*sqrt(1 + rho_R^2)
}

##### Initialize arrays #####

# Set historical FM to 0
Fb <- 0

# Fishing effort stays constant
E2 <- array(rep(1, eq_time), c(1, eq_time, 1, 1))

# Initialize FM and depletion levels
FM_values <- seq(from = 0, to = 2, by = 0.01)
fn <- length(FM_values)

# Initialize depletion vector
dep <- rep(0, fn)

# Initialize population size and catch arrays
# Dimensions = age * 1 * time * 1
N2 <- catch2 <- array(rep(0, n*eq_time), c(n, 1, eq_time, 1, 1))

# Initialize biomass, SSB, and recruitment error
# Dimensions = 1 * time * 1
SSB2 <- biomass2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
abundance_all2 <- abundance_mature2 <- array(rep(0, eq_time), 
                                             c(1, eq_time, 1, 1))

rec_biomass <- array(rep(NA, fn*eq_time), c(fn, eq_time))

SAD <- equilibrium_SAD(rec_age, max_age, n, W, R0, Mat, h, B0, sigma_R, Fb, 
                       S, M, eq_time = 150, m, stochasticity = F, rho_R, 
                       nat_mortality = M, recruitment_mode, A)

# Enter N, abundance, catch, and biomasses for time = 1 to rec_age
for (t in 1:rec_age) {
  N2[, 1, t, 1, 1] <- SAD
  biomass2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1] * W)
  catch2[, 1, t, 1, 1] <- rep(0, n)
  SSB2[1, t, 1, 1] <- sum(N2[, 1, t - rec_age, 1, 1]*W*Mat)
  abundance_all2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1])
  abundance_mature2[1, t, 1, 1] <- sum(N2[m:max_age - 1, 1, t, 1, 1])
} 

# Initialize FM matrix
FM2 <- array(rep(0, n*eq_time), c(n, 1, eq_time, 1, 1))

# Step population forward in time with set fishing level
for (t in (rec_age + 1):eq_time) {
  
  PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, rec_age, max_age, n, SSB2, 
                     N2, W, Mat, A = 1, R0, h, B0, Eps2, sigma_R, Fb = 0, E2, 
                     S, NM, FM2, m, abundance_all2, abundance_mature2, 
                     biomass2, fishing = T, nat_mortality = M, recruitment_mode)
  
  FM2[, 1, t, 1, 1]               <- rep(0, n)
  N2[, 1, t, 1, 1]                <- PD[[2]]
  abundance_all2[1, t, 1, 1]      <- PD[[3]]
  abundance_mature2[1, t, 1, 1]   <- PD[[4]]
  biomass2[1, t, 1, 1]            <- PD[[5]]
  SSB2[1, t, 1, 1]                <- PD[[6]]
  
}

# Calculate final biomass given zero fishing
FM0_biomass <- biomass2[1, eq_time, 1, 1]

# Substitute in values for Fb to get depletion level
for (i in 2:fn) { 
  
  FM2 <- array(rep(FM_values[i], n*eq_time), c(n, 1, eq_time, 1, 1))
  
  # Step population forward in time with set fishing level
  for (t in (rec_age + 1):eq_time) {
    
    PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, rec_age, max_age, n, SSB2, 
                       N2, W, Mat, A = 1, R0, h, B0, Eps2, sigma_R, Fb = 0, E2, 
                       S, NM, FM2, m, abundance_all2, abundance_mature2, 
                       biomass2, fishing = T, nat_mortality = M, 
                       recruitment_mode)
    
    FM2[, 1, t, 1, 1]               <- rep(FM_values[i], n)
    N2[, 1, t, 1, 1]                <- PD[[2]]
    abundance_all2[1, t, 1, 1]      <- PD[[3]]
    abundance_mature2[1, t, 1, 1]   <- PD[[4]]
    biomass2[1, t, 1, 1]            <- PD[[5]]
    SSB2[1, t, 1, 1]                <- PD[[6]]
    
  }
  
  dep[i] <- 1 - (biomass2[1, eq_time, 1, 1] / FM0_biomass)

}

closest_FM <- FM_values[which.min(abs(dep - true_dep))] 

plot(FM_values, dep)
abline(v = closest_FM, col = 'red')
abline(h = true_dep, col = 'green')
print(closest_FM)
