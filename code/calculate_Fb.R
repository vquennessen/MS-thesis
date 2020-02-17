setwd('~/Projects/MS-thesis/code')

devtools::install_git('https://github.com/vquennessen/densityratio.git')

library(densityratio)

# # black rockfish (OR) 2015
# historical_FM(Species = 'BR_OR_2015')
# 
# # cabezon (OR) 2019
# historical_FM(Species = 'CAB_OR_2019')
# 
# # lingcod (OR and WA) 2017
# historical_FM(Species = 'LING_OW_2017')
# 
# # canary rockfish (OR) 2015
# historical_FM(Species = 'CR_OR_2015')
# 
# # english sole (Northern CA, OR, and WA) 2013
# historical_FM(Species = 'ES_COW_2013')

species_list <- c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015', 
                  'ES_COW_2013')

eq_time <- 150
R0 <- 1e+5
Stochasticity <- FALSE
Recruitment_mode <- 'pool'

for (i in 1:length(species_list)) {
  
  par <- parameters(species_list[i])
  
  Max_age                <- par[[1]]        # maximum age
  M                      <- par[[2]]        # natural mortality
  Rec_age                <- par[[3]]        # age at recruitment
  WA  <- par[[4]];  WB   <- par[[5]]        # weight at length parameters (f)
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
  SM                     <- par[[24]]       # is selectivity = maturity at age?
  Fleets                 <- par[[25]]       # fishery fleet names
  Alpha                  <- par[[26]]       # slope for upcurve
  Beta                   <- par[[27]]       # slope for downcurve
  F_fin                  <- par[[28]]       # F_fin for fishery, 0 if asymptotic
  A50_up                 <- par[[29]]       # L50 for upcurve
  A50_down               <- par[[30]]       # L50 for downcurve
  Cf                     <- par[[31]]       # fraction of fishery caught / fleet
  
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
  B0 <- R0/Phi
  # selectivity at age
  if (SM == FALSE) {
    S <- selectivity(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets, A50_up,
                     A50_down, Alpha, F_fin, Beta, Cf)
  } else { S <- Mat }
  
  plot(ages, S - Mat, main = species_list[i], ylim = c(-1, 1))
  
}
