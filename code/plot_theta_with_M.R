# plot total theta for static vs. transient DRs over time

# load any necessary libraries
library(ggplot2)
library(patchwork)
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
num_sims <- 2
data_folder <- 'None'
figures_folder <- 'None/relative'
###############################################################################

# species to compare
species_list <- c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')
Names <- c('Black Rockfish', 'Cabezon', 'Lingcod', 'Canary Rockfish')

# set variables
A = 5
MPA = 3
Time2 = 20
Final_DRs <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules <- 1:6
R0 <- 1e5

# dimensions
PD = 0.25
Error = 0.05
estimates <- c('Low', 'True', 'High')
ENM = 2
sample_size = num_sims

nC <- length(Control_rules)
nT <- Time2 + 1
nE <- length(estimates)
nS <- length(species_list)
nF <- length(Final_DRs)

theta_df <- data.frame(Species = rep(Names, each = nE*2*nF*nT),
                       Type = rep(c('Static', 'Transient'), 
                                  each = nE*nF*nT, times = nS),
                       Estimate = rep(estimates, times = nS*2, each = nF*nT), 
                       FDR = rep(Final_DRs, each = nT, times = 2*nS), 
                       Time = rep(0:Time2, times = nF*2*nS),
                       Theta = rep(0, nS*2*nF*nT))

for (s in 1:length(species_list)) {
  
  # load objects
  load(paste('~/Projects/MS-thesis/data/', data_folder, '/', 
             species_list[s], '/', num_sims, '_N.Rda', sep = ''))
  
  ##### relative biomass and median, upper, and lower limits  #####
  
  # calculate stable age distribution
  par <- parameters(species_list[s])
  
  Max_age                <- par[[1]]        # maximum age
  M                      <- par[[2]]        # natural mortality
  Rec_age                <- par[[3]]        # age at recruitment
  WA <- par[[4]];  WB    <- par[[5]]        # weight at length parameters (f)
  A1 <- par[[6]];  L1    <- par[[7]]        # growth parameters (f)
  A2 <- par[[8]];  L2    <- par[[9]]
  K                      <- par[[10]]
  L50                    <- par[[11]]       # length at 50% maturity
  K_mat                  <- par[[12]]       # slope of maturity curve
  H                      <- par[[13]]       # steepness
  Phi                    <- par[[14]]       # unfished recruits per spawner
  Sigma_R                <- par[[15]]       # recruitment standard deviation
  Rho_R                  <- par[[16]]       # recruitment autocorrelation
  AMP                    <- par[[17]]       # adult movement proportion
  D                      <- par[[18]]       # depletion
  Fb                     <- par[[19]]       # fishing mortality to cause D
  Fleets                 <- par[[23]]       # fishery fleet names
  Alpha                  <- par[[24]]       # slope for upcurve
  Beta                   <- par[[25]]       # slope for downcurve
  F_fin                  <- par[[26]]       # F_fin for fishery, 0 if asymptotic
  A50_up                 <- par[[27]]       # A50 for upcurve
  A50_down               <- par[[28]]       # A50 for downcurve
  Cf                     <- par[[29]]       # fraction of fishery caught / fleet
  
  ages <- Rec_age:Max_age
  n <- length(ages)
  L <- length_age(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = FALSE)
  W <- weight(L, WA, WB)
  Mat <- maturity(Rec_age, Max_age, K_mat, L, L50)
  S <- selectivity(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets, A50_up,
                   A50_down, Alpha, F_fin, Beta, Cf)
  B0 <- R0 / Phi
  A50_mat <- ages[min(which(Mat > 0.5))]
  SAD <- stable_AD(Rec_age, Max_age, W, R0, Mat, H, B0, Sigma_R, Fb, S, M, 
                   eq_time = 150, A50_mat, Stochasticity = FALSE, Rho_R, 
                   Recruitment_mode = 'pool', A)
  
  # pull out sample sims and calculate median across num_sims
  A_sample <- array(rep(0, n*nT*nC*nF*num_sims), c(n, nT, nC, nF, num_sims))
  A_medians <- array(rep(0, n*nT*nC*nF), c(n, nT, nC, nF))
  
  for (i in 1:n) {
    A_sample[i, , , , ] <- colSums(sims_N[i, , , , , ])
    for (t in 1:nT) {
      for (cr in 1:nC) {
        for (fdr in 1:nF) {
          A_medians[i, t, cr, fdr] <- median(A_sample[i, t, cr, fdr, ])
        }
      }
    }
  }
  
  for (type in 1:2) {
    for (e in 1:nE) {
      for (fdr in 1:nF) {
        for (t in 1:nT) {
          num <- SAD %*% A_medians[, t, (type - 1)*3 + e, fdr]
          norm1 <- norm(matrix(SAD), type = 'F')
          norm2 <- norm(matrix(A_medians[, t, (type - 1)*3 + e, fdr]), 
                        type = 'F')
          denom <- norm1 * norm2
          i <- (s - 1)*2*nE*nF*nT + (type - 1)*nE*nF*nT + (e - 1)*nF*nT + 
            (fdr - 1)*nT + t
          theta_df$Theta[i] <- as.numeric(acos(num / denom))
        }
      }
    }
  }
  
}

nI <- 2*nE*nF*nT

dfA <- theta_df[1:nI, ]
dfB <- theta_df[(nI + 1):(2*nI), ]
dfC <- theta_df[(2*nI + 1):(3*nI), ]
dfD <- theta_df[(3*nI + 1):(4*nI), ]

dfA <- subset(dfA, FDR < 0.7)
dfB <- subset(dfB, FDR < 0.7)
dfC <- subset(dfC, FDR < 0.7)
dfD <- subset(dfD, FDR < 0.7)


# plotting parameters
y1 <- 0
y2 <- 0.7
jitter_height <- 0.005
size1 <- 1
size2 <- 0.5

# plot it
A <- ggplot(data = dfA, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = Type)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle('Black Rockfish') +
  ylab('Theta') +
  ylim(y1, y2) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

B <- ggplot(data = dfB, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = Type)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle('Cabezon') +
  ylim(y1, y2) +
  theme(axis.text = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = 'none') 

C <- ggplot(data = dfC, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = Type)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle('Lingcod') +
  xlab('Years since reserve implemented') +
  ylab('Theta') +
  ylim(y1, y2) +
  theme(legend.position = 'none')

D <- ggplot(data = dfD, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = Type)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle('Copper Rockfish') +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ylim(y1, y2) +
  xlab('Years since reserve implemented') +
  theme(plot.margin = unit(c(0, 80, 0, 0), 'pt')) +
  labs(linetype = 'Type', color = 'FDR') +
  theme(legend.position = c(1.25, 1))

##### patch all the figures together #####
patch2 <- (A + B) / (C + D)
thing2 <- patch2 + plot_annotation(
  title = 'Distance of total population from stable age distribution')

ggsave(thing2, filename = 'M_relative_theta.png',
       path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                    figures_folder, sep = ''))
