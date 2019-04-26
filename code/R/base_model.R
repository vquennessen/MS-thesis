#' Runs age-structured population dynamics biological sub-model. 

##### Load life history characteristics for species ############################

### Black rockfish (Sebastes melanops) ###

max_age <- 40                           ## maximum age
M <- 0.17                               ## natural mortality
rec_age <- 3                            ## age at recruitment
af <- 2.6e-5; bf <- 2.88                ## weight at length parameters (f)
am <- 2.58e-5; bm <- 2.89               ## weight at length parameters (m)
a1f <- 1; L1f <- 20.32;                 ## growth parameters (f)
a2f <- 40; L2f <- 49.67; Kf <- 0.21  
a1m <- 1; L1m <- 17.47;                 ## growth parameters (m)
a2m <- 40; L2m <- 43.27; Km <- 0.34  
L50 <- 43.69                            ## length at 50% maturity
k_mat <- -0.66                          ## slope of maturity curve
ldp <- 0.1                              # larval drift proportion
R0 <- 3666                              ## unfished recruitment
h <- 0.77                               ## steepness
phi <- 1.1                              # unfished recruits per spawner
sigma_R <- 0.5                          ## recruitment standard deviation
rho_R <- 0                              # recruitment autocorrelation
p <- 0.1                                # adult movement proportion
D <- 0.604                              ## depletion
Fb <- 0.2                               # fishing mortality to cause D
r <- 0.77                               # Proportion of positive transects 
                                        #       in PISCO monitoring data
x <- 15.42                              # mean of positive transects
sp <- 16.97                             # std of positive transects
fleet_alpha <- 1                        # selectivity paramter
fleet_beta <- 1                         # selectivity paramter
# TODO: get actual alpha and beta values / fleet values???


##### Population Dynamics ######################################################

s <- 5            # number of areas
t <- 50           # number of timesteps (years)
E <- 0.10         # nominal fishing effort in each area 
# TODO: update E value

# Initialize arrays 
n <- rec_age:max_age
N <- array(rep(0, length(n)*s*t), c(length(n), s, t))
B <- array(rep(0, s*t), c(s, t))

# Length at age
# Based on Babcock & MacCall (2011): Eq. (10)
age <- array(rep(0, max_age), max_age)
L_inf <- L1f + (L2f - L1f)/(1 - exp(-1*Kf*(a2f - a1f)))
L <- L_inf + (L1f - L_inf)*exp(-1*Kf*(age - a1f))

# Weight at age
# Based on Babcock & MacCall (2011): Eq. (11)
W <- af*L^bf

# Maturity at length
M <- array(rep(0, max_age), max_age)
M <- (1)/(1 + exp(k_mat*(L - L50)))

# Based on Babcock & MacCall (2011): Eq. (4)
epsilon <- array(rep(0, t), t)
nu <- rnorm(t, 0, sigma_R)
for (i in 2:t) {
  epsilon[i] <- rho_R*epsilon[i-1] + nu[i]*sqrt(1 + rho_R^2)
}

# Spawning Stock Biomass


# Recruitment
# Based on Babcock & MacCall (2011): Eq. (3)
R <- array(rep(0, s*t), c(s, t))

# Catchability
# Based on Babcock & MacCall (2011): Eq. (6)
q <- (s*Fb)/(s*E)

# Selectivity
# Based on Babcock & MacCall (2011): Eq. (8)
Sp1 <- (1)/(1 + exp(-1*fleet_alpha[1]*(L - L50)))
# Based on Babcock & MacCall (2011): Eq. (9)
Sp2 <- 1 - (1 - Ffin)/(1 + exp(-1*fleet_beta[1]*(L - L50)))
# TODO: Figure out what Ffin is
# TODO: Figure out selectivity, alpha and beta values

# Selectivity at Age
# Based on Babcock & MacCall (2011): Eq. (7)
Sp <- array(rep(0, n*length(fleet_alpha)), c(n, length(fleet_alpha)))

# Fishing Mortality
# Based on Babcock & MacCall (2011): Eq. (5)
f <- array(rep(0, length(n)*s*t), c(length(n), s, t))

# Population size for different ages
# Based on Babcock & MacCall (2011): Eq. (1)

