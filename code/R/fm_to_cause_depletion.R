# set working directory
setwd("~/Projects/DensityRatio/code/R")

# clear environment
rm(list = ls())

##### specify parameters for this run #####
Species <- 'BR2003'
eq_time <- 150
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

Max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
Rec_age                <- par[[3]]        # age at recruitment
a   <- par[[4]];  b    <- par[[5]]        # weight at length parameters (f)
A1  <- par[[6]];  L1   <- par[[7]]        # growth parameters (f)
A2  <- par[[8]];  L2   <- par[[9]] 
K   <- par[[10]]  
L50                    <- par[[11]]       # length at 50% maturity
K_mat                  <- par[[12]]       # slope of maturity curve
H                      <- par[[14]]       # steepness
Phi                    <- par[[15]]       # unfished recruits per spawner
Sigma_R                <- par[[16]]       # recruitment standard deviation
Rho_R                  <- par[[17]]       # recruitment autocorrelation
                                          #       in PISCO monitoring data
D                      <- par[[19]]
SP                     <- par[[23]]       # std of positive transects

####### selectivity parameters #######
Fleets                 <- par[[26]]       # fishery fleet names
Alpha                  <- par[[27]]       # slope for upcurve
Beta                   <- par[[28]]       # slope for downcurve
F_fin                  <- par[[29]]       # F_fin for fishery, 0 if asymptotic
A50_up                 <- par[[30]]       # L50 for upcurve
A50_down               <- par[[31]]       # L50 for downcurve
Cf                     <- par[[32]]       # fraction of fishery caught / fleet

##### Calculate set values #####
Age <- rec_age:max_age                            # applicable ages
Num <- length(age)                                # number of age bins
L   <- length_at_age(Rec_age, Max_age, A1, L1, 
                     A2, L2, K, All_ages = F)     # length at age
W <- weight_at_age(L, WA, WB)                     # weight at age
Mat <- fraction_mature_at_age(Num, K_mat, L, L50) # maturity at age
m <- age[min(which(Mat > 0.5))]                   # age at 50% mature
B0 <- R0/phi                                      # unfished SSB
S <- selectivity_at_age(Fleets, L, Max_age, Rec_age, Alpha, A50_up, 
                        A50_down, F_fin, Beta, Num, Cf, Age)

# Recruitment normal variable
# Dimensions = area * timeT * CR
if (Stochasticity == T) {
  nuR2 <- array(rnorm(eq_time, 0, Sigma_R), c(1, eq_time, 1, 1))
} else if (stochasticity == F) {
  nuR2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
}

# initialize epsilon vector
Eps2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
# eps[1]
Eps2[1, 1, 1, 1] <- nuR2[1, 1, 1, 1]*sqrt(1 + Rho_R^2)
# fill in rest of epsilon vector
for (t in 2:eq_time) {
  Eps2[1, t, 1, 1] <- Rho_R*Eps2[1, t-1, 1, 1] + 
    nuR2[1, t, 1, 1]*sqrt(1 + Rho_R^2)
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

SAD <- equilibrium_SAD(Rec_age, Max_age, Num, W, R0, Mat, H, B0, Sigma_R, Fb, 
                       S, M, eq_time = 150, A50_mat, Stochasticity = F, Rho_R, 
                       Nat_mortality = M, Recruitment_mode, A)

# Enter N, abundance, catch, and biomasses for time = 1 to rec_age
for (t in 1:rec_age) {
  N2[, 1, t, 1, 1] <- SAD
  biomass2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1] * W)
  catch2[, 1, t, 1, 1] <- rep(0, n)
  SSB2[1, t, 1, 1] <- sum(N2[, 1, t - Rec_age, 1, 1]*W*Mat)
  abundance_all2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1])
  abundance_mature2[1, t, 1, 1] <- sum(N2[A50_mat:Max_age - 1, 1, t, 1, 1])
} 

# Initialize FM matrix
FM2 <- array(rep(0, Num*eq_time), c(Num, 1, eq_time, 1, 1))

# Step population forward in time with set fishing level
for (t in (Rec_age + 1):eq_time) {
  
  PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, Rec_age, Max_age, Num, SSB2, 
                     N2, W, Mat, A = 1, R0, H, B0, Eps2, Sigma_R, Fb = 0, E2, 
                     S, NM, FM2, A50_mat, abundance_all2, abundance_mature2, 
                     biomass2, Fishing = T, Nat_mortality = M, Recruitment_mode)
  
  FM2[, 1, t, 1, 1]               <- rep(0, Num)
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
  
  FM2 <- array(rep(FM_values[i], Num*eq_time), c(Num, 1, eq_time, 1, 1))
  
  # Step population forward in time with set fishing level
  for (t in (Rec_age + 1):eq_time) {
    
    PD <- pop_dynamics(a = 1, t, cr = 1, nm = 1, Rec_age, Max_age, Num, SSB2, 
                       N2, W, Mat, A = 1, R0, H, B0, Eps2, Sigma_R, Fb = 0, E2, 
                       S, NM, FM2, A50_mat, abundance_all2, abundance_mature2, 
                       biomass2, Fishing = T, Nat_mortality = M, 
                       Recruitment_mode)
    
    FM2[, 1, t, 1, 1]               <- rep(FM_values[i], Num)
    N2[, 1, t, 1, 1]                <- PD[[2]]
    abundance_all2[1, t, 1, 1]      <- PD[[3]]
    abundance_mature2[1, t, 1, 1]   <- PD[[4]]
    biomass2[1, t, 1, 1]            <- PD[[5]]
    SSB2[1, t, 1, 1]                <- PD[[6]]
    
  }
  
  dep[i] <- 1 - (biomass2[1, eq_time, 1, 1] / FM0_biomass)

}

closest_FM <- FM_values[which.min(abs(dep - D))] 

plot(FM_values, dep)
abline(v = closest_FM, col = 'red')
abline(h = D, col = 'green')
print(closest_FM)
