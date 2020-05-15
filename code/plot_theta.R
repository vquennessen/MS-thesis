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
folder <- 'None/relative'
###############################################################################

# species to compare
species_list <- c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')
Names <- c('Black Rockfish', 'Cabezon', 'Lingcod', 'Canary Rockfish')

# sample from simulations
indices <- sample(1:num_sims, sample_size, replace = FALSE)

# set variables
A = 5
MPA = 3
R0 = 1e5
Time2 = 20
Final_DRs = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:2)

# dimensions
sample_size = num_sims
PD = 0.25

nT <- Time2 + 1
nS <- length(species_list)
nF <- length(Final_DRs)
nA <- MPA
nC <- length(Control_rules)

theta_df <- data.frame(Species = rep(species_list, each = nC*nF*nT),
                       Type = rep(c('Static', 'Transient'), 
                                  each = nF*nT, times = nS),
                       FDR = rep(Final_DRs, each = nT, times = nC*nS), 
                       Time = rep(0:Time2, times = nF*nC*nS),
                       Theta = rep(0, nS*nC*nF*nT))

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
  
  
  # initialize theta values
  theta <- array(rep(NA, nT*nC*nF), c(nT, nC, nF))
  
  for (cr in 1:nC) {
    for (fdr in 1:nF) {
      for (t in 1:nT) {
        num <- SAD %*% A_medians[, t, cr, fdr]
        norm1 <- norm(matrix(SAD), type = 'F')
        norm2 <- norm(matrix(A_medians[, t, cr, fdr]), type = 'F')
        denom <- norm1 * norm2
        type <- ifelse(cr < 4, 1, 2)
        i <- (s - 1)*nC*nF*nT + (cr - 1)*nF*nT + (fdr - 1)*nT + t
        theta_df$Theta[i] <- as.numeric(acos(num / denom))
      }
    }
  }
  
}

nI <- nC*nF*nT

dfA <- theta_df[1:nI, ]
dfB <- theta_df[(nI + 1):(2*nI), ]
dfC <- theta_df[(2*nI + 1):(3*nI), ]
dfD <- theta_df[(3*nI + 1):(4*nI), ]

# plotting parameters
y1 <- 0
y2 <- 0.7
jitter_height <- 0.005
ltype1 <- 2
ltype2 <- 1

# plot it
A <- ggplot(data = dfA, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_linetype_manual(values = c(ltype1, ltype2)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle('Black Rockfish') +
  ylab('Theta') +
  ylim(y1, y2) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

B <- ggplot(data = dfB, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_linetype_manual(values = c(ltype1, ltype2)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle('Cabezon') +
  ylim(y1, y2) +
  theme(axis.text = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = 'none') 

C <- ggplot(data = dfC, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_linetype_manual(values = c(ltype1, ltype2)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle('Lingcod') +
  xlab('Years since reserve implemented') +
  ylab('Theta') +
  ylim(y1, y2) +
  theme(legend.position = 'none')

D <- ggplot(data = dfD, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_linetype_manual(values = c(ltype1, ltype2)) +
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

ggsave(thing2, filename = 'all_species_theta.png',
       path = paste('C:/Users/Vic/Google Drive/OSU/Thesis/figures/', folder, sep = ''))
