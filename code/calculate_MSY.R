#' Calculate precautionary MSY for species

# Black Rockfish MSY = 40% Unfished Biomass
# Cabezon MSY = 40% Unfished Biomass
# Lingcod MSY = 40% Unfished Spawning Stock Biomass
# Canary Rockfish MSY = 40% Unfished Biomass

# load libraries
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

##### model parameters #########################################################

Species <- 'CR_OR_2015'
# Species <- 'BR_OR_2015'
# Species <- 'LING_OW_2017'
# Species <- 'CAB_OR_2019'

eq_time <- 250


##### stable age distribution ##################################################



##### load species parameters #####
par <- parameters(Species)

Max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
Rec_age                <- par[[3]]        # age at recruitment
WA  <- par[[4]];  WB   <- par[[5]]        # weight at length parameters (f)
A1  <- par[[6]];  L1   <- par[[7]]        # growth parameters (f)
A2  <- par[[8]];  L2   <- par[[9]]
K   <- par[[10]]
L50                    <- par[[11]]       # length at 50% maturity
K_mat                  <- par[[12]]       # slope of maturity curve
H                      <- par[[13]]       # steepness
Sigma_R                <- par[[14]]       # recruitment standard deviation
Rho_R                  <- par[[15]]       # recruitment autocorrelation
#       in PISCO monitoring data
D                      <- par[[17]]
SP                     <- par[[21]]       # std of positive transects
Fleets                 <- par[[22]]       # fishery fleet names
Alpha                  <- par[[23]]       # slope for upcurve
Beta                   <- par[[24]]       # slope for downcurve
F_fin                  <- par[[25]]       # F_fin for fishery, 0 if asymptotic
A50_up                 <- par[[26]]       # L50 for upcurve
A50_down               <- par[[27]]       # L50 for downcurve
Cf                     <- par[[28]]       # fraction of fishery caught / fleet

##### Calculate set values #####
ages <- Rec_age:Max_age                            # applicable ages
num <- length(ages)                                # number of age bins
# length at age
L   <- length_age(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = F)
# weight at age
W <- weight(L, WA, WB)
# fraction mature at age
Mat <- maturity(Rec_age, Max_age, K_mat, L, L50)
# age at 50% mature
A50_mat <- ages[min(which(Mat > 0.5))]
# unfished biomass
R0 <- 1000
B0 <- R0*W[1]
# selectivity at age
S <- selectivity(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets, A50_up,
                 A50_down, Alpha, F_fin, Beta, Cf)

# # Recruitment error = 0 without Recruitment_Var
# Eps2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))

# stable age distribution
SAD <- stable_AD(Rec_age, Max_age, W, R0, Mat, H, B0, Sigma_R, Fb = 0, S, M, 
                 eq_time, A50_mat, Recruitment_Var = FALSE, Rho_R, 
                 Recruitment_mode = 'pool', A = 5)

##### Initialize arrays ########################################################

# Fishing effort stays constant
E2 <- array(rep(1, eq_time), c(1, eq_time, 1, 1))

# Initialize FM and depletion levels
FM_values <- seq(from = 0, to = 1, by = 0.01)
fn <- length(FM_values)

# Initialize depletion vectors
dep_SSB <- rep(0, fn)
dep_B <- rep(0, fn)

# Initialize population size and catch arrays
# Dimensions = age * 1 * time * 1
N2 <- catch2 <- array(rep(0, num*eq_time), c(num, 1, eq_time, 1,  1))

# Initialize biomass, SSB, and abundance
# Dimensions = 1 * time * 1
SSB2 <- biomass2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1))
abundance2 <- array(rep(0, eq_time), c(1, eq_time, 1, 1, 1))

# Initialize catch-at-age matrix
# Dimensions = age * area * time * CR * FDR values (3)
Catch <- array(rep(0, num*1*TimeT*1*1), c(num, 1, TimeT, 1, 1))

# Initialize yield matrix
# Dimensions = area * time * CR * FDR values (3)
Yield <- array(rep(0, 1*TimeT*1*1), c(1, TimeT, 1, 1))

# initial spawning stock biomass / biomass with no fishing
FM0_SSB <- sum(W*SAD*Mat)
FM0_B <- sum(W*SAD)

for (t in 1:Rec_age) {
  N2[, 1, t, 1, 1] <- SAD
  abundance2[1, t, 1, 1, 1] <- sum(N2[, 1, t, 1, 1])
  biomass2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1] * W)
  SSB2[1, t, 1, 1] <- sum(N2[, 1, t, 1, 1]*W*Mat)
}

# Substitute in values for Fb to get depletion level
for (i in 2:fn) {
  
  FM2 <- array(rep(FM_values[i], num*eq_time), c(num, 1, eq_time, 1, 1))
  
  # Step population forward in time with set fishing level
  for (t in (Rec_age + 1):eq_time) {
    
    # recruitment
    R <- recruitment(t, cr = 1, fdr = 1, SSB2, A = 1, R0, H, B0, Eps2,
                     Sigma_R, Rec_age, Recruitment_mode, LDP)
    
    PD <- pop_dynamics(t, cr = 1, fdr = 1, Rec_age, Max_age, SSB2, N2, W, Mat,
                       A = 1, Fb = 0, E2, S, FM2, A50_mat, biomass2,
                       abundance2, Fishing = T, Nat_mortality = M, R)
    
    FM2[, , t, 1, 1]               <- rep(FM_values[i], num)
    N2[, , t, 1, 1]                <- PD[[2]]
    biomass2[, t, 1, 1]            <- PD[[3]]
    SSB2[, t, 1, 1]                <- PD[[4]]
    abundance2[, t, 1, 1, 1]       <- PD[[5]]
    
    # fishing
    Catch[, 1, t, 1, 1] <- catch(t, cr, fdr, FM, Nat_mortality, N, A = 1, Fb, 
                                 E, Catch)
    Yield[, t, 1, 1] <- colSums(Catch[, 1, t, cr, fdr]*W)
    
    }
  
  dep_SSB[i] <- SSB2[1, eq_time, 1, 1] / FM0_SSB
  dep_B[i] <- B2[1, eq_time, 1, 1] / FM0_B
  
}

FMSY_SSB <- FM_values[which.min(abs(dep_SSB - 0.4))]
FMSY_B <- FM_values[which.min(abs(dep_B - 0.4))]

# MSY by SSB
plot(FM_values, dep_SSB, main = Species, ylim = c(0, 1),
     ylab = 'Depletion', xlab = 'FM value')
abline(v = FMSY_SSB, col = 'red')
abline(h = 0.4, col = 'green')

# MSY by B
plot(FM_values, dep_B, main = Species, ylim = c(0, 1),
     ylab = 'Depletion', xlab = 'FM value')
abline(v = FMSY_B, col = 'red')
abline(h = 0.4, col = 'green')
