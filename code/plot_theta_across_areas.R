# plot theta for static vs. transient DRs over time

# load any necessary libraries
library(ggplot2)
library(patchwork)
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

# species to compare
species_list <- c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')

# set variables
num_sims <- 2
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
Control_rules = c(1:6)

# dimensions
TimeT <- Time1 + Time2
sample_size = num_sims
PD = 0.25
plot_individual_runs = FALSE
Error = 0.05
estimates <- c('Low', 'True', 'High')
ENM = 2
R0 = 1e5

nT <- Time2 + 1
nE <- length(estimates)
nS <- length(species_list)
nF <- length(Final_DRs)
nA <- MPA

base_df <- data.frame(Species = rep(species_list, each = 2*nE*nF*nT),
                      Type = rep(c('Static', 'Transient'), 
                                         each = nE*nF*nT, times = nS),
                      Estimate = rep(estimates, each = nT*nF, times = 2*nS), 
                      FDR = rep(Final_DRs, each = nT, times = nE*2*nS), 
                      Time = rep(0:Time2, times = nF*nE*2*nS),
                      Theta = rep(0, nS*2*nE*nF*nT))

for (s in 1:length(species_list)) {
  
  # load objects
  load(paste('~/Projects/MS-thesis/data/no_stochasticity/all_FDR_values/', 
             species_list[s], '/', num_sims, '_N.Rda', sep = ''))
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
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
  
  # nat_mortality
  Nat_mortality <- c(M - Error, M, M + Error)
  NM <- length(Nat_mortality)
  nC <- length(Control_rules)
  
  # pull out sample sims
  A_sample <- sims_N[, , , , , indices]
  
  # initialize theta values
  theta <- array(rep(NA, nA*nT*nC*nF*num_sims), c(nA, nT, nC, nF, num_sims))
  
  for (a in 1:nA) {
    for (t in 1:nT) {
      for (cr in 1:nC) {
        for (fdr in 1:nF) {
          for (ns in 1:num_sims) {
            num <- SAD %*% sims_N[, a, t, cr, fdr, ns]
            norm1 <- norm(matrix(SAD), type = 'F')
            norm2 <- norm(matrix(sims_N[, a, t, cr, fdr, ns]), type = 'F')
            denom <- norm1 * norm2
            theta[a, t, cr, fdr, ns] <- acos(num / denom)
          }
        }
      }
    }
  }

##### far - initialize and fill in data frame #####

far_df <- base_df
near_df <- base_df
inside_df <- base_df
  
for (c in 1:2) {
  for (e in 1:nE) {
    for (fdr in 1:nF) {
      for (t in 1:nT) {
        i <- (s-1)*2*nE*nF*nT + (c-1)*nE*nF*nT + (e-1)*nF*nT + (fdr-1)*nT + t
        far_df$Theta[i] <- median(theta[1, t, 3*(c - 1) + e, fdr, ])
        near_df$Theta[i] <- median(theta[2, t, 3*(c - 1) + e, fdr, ])
        inside_df$Theta[i] <- median(theta[3, t, 3*(c - 1) + e, fdr, ])
      }
    }
  }
}

# plotting parameters
y1 = 0
y1.5 = 0
y2 = 0.8
y2.5 = 0.31
jitter_height1 <- 0.01
jitter_height2 <- 0.001

sp_far <- far_df[which(far_df$Species == species_list[s]), ]

# legend below all three graphs side by side
# plot far from reserve
# FAR <- ggplot(data = sp_far, aes(x = Time, y = Theta, 
#                                   color = as.factor(FDR), 
#                                   linetype = as.factor(Estimate), 
#                                   size = as.factor(Type))) +
#   geom_line(position = position_jitter(w = 0, h = 0.002)) +
#   scale_size_manual(values = c(0.5, 1.5), guide = 'none') +
#   ylim(y1, y2) +
#   ggtitle('Far from Reserve') +
#   theme(axis.title.x = element_blank()) +
#   labs(color = 'FDR', linetype = 'M', size = 'Type') +
#   theme(legend.position = 'bottom')  +
#   scale_colour_hue(guide = 'none')

FAR <- ggplot(data = sp_far, aes(x = Time, y = Theta, 
                                 color = as.factor(FDR), 
                                 linetype = as.factor(Estimate), 
                                 size = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height1)) +
  scale_size_manual(values = c(0.5, 1.5)) +
  ylim(y1, y2) +
  ggtitle('Far from Reserve') +
  theme(axis.title.x = element_blank()) +
  labs(color = 'FDR', linetype = 'M', size = 'Type') +
  theme(legend.position = 'none')
# near DFs
sp_near <- near_df[which(near_df$Species == species_list[s]), ]

# legend below all three graphs side by side
# plot near reserve
# NEAR <- ggplot(data = sp_near, aes(x = Time, y = Theta, 
#                                 color = as.factor(FDR), 
#                                 linetype = as.factor(Estimate), 
#                                 size = as.factor(Type))) +
#   geom_line(position = position_jitter(w = 0, h = 0.002)) +
#   
#   scale_size_manual(values = c(0.5, 1.5), guide = 'none') +
#   ylim(y1, y2) +
#   ggtitle('Near Reserve') +
#   xlab('Years since reserve implementation') +
#   theme(axis.title.y = element_blank()) +
#   theme(axis.text.y = element_blank()) +
#   labs(color = 'FDR', linetype = 'M', size = 'Type') +
#   theme(legend.position = 'bottom')  +
#   scale_linetype(guide = 'none')

NEAR <- ggplot(data = sp_near, aes(x = Time, y = Theta, 
                                   color = as.factor(FDR), 
                                   linetype = as.factor(Estimate), 
                                   size = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height1)) +
  
  scale_size_manual(values = c(0.5, 1.5)) +
  ylim(y1, y2) +
  ggtitle('Near Reserve') +
  xlab('Years since reserve implementation') +
  theme(axis.title = element_blank()) +
  theme(axis.text.y = element_blank()) +
  labs(color = 'FDR', linetype = 'M', size = 'Type') +
  theme(legend.position = 'none')
# plot inside reserve
sp_inside <- inside_df[which(inside_df$Species == species_list[s]), ]

# legend below all three graphs side by side
# INSIDE <- ggplot(data = sp_inside, aes(x = Time, y = Theta, 
#                                      color = as.factor(FDR), 
#                                      linetype = as.factor(Estimate), 
#                                      size = as.factor(Type))) +
#   geom_line(position = position_jitter(w = 0, h = 0.0002)) +
#   ylim(y1.5, y2) +
#   scale_size_manual(values = c(0.5, 1.5)) +
#   ggtitle('Inside Reserve') +
#   theme(axis.title.x = element_blank()) +
#   theme(axis.title.y = element_blank()) +
#   labs(color = 'FDR', linetype = 'M', size = 'Type') +
#   theme(legend.position = 'bottom') +
#   scale_colour_hue(guide = 'none') + scale_linetype(guide = 'none')

INSIDE <- ggplot(data = sp_inside, aes(x = Time, y = Theta, 
                                       color = as.factor(FDR), 
                                       linetype = as.factor(Estimate), 
                                       size = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height2)) +
  ylim(y1.5, y2.5) +
  scale_size_manual(values = c(0.5, 1.5)) +
  ggtitle('Inside Reserve') +
  xlab('Time (years) since reserve implementation') +
  labs(color = 'FDR', linetype = 'M', size = 'Type') +
  theme(legend.position = c(1.5, 0.4), legend.box = 'horizontal')

##### patch all the figures together #####
patch <- (FAR + NEAR) / (INSIDE + plot_spacer())
# patch <- FAR + NEAR + INSIDE
thing <- patch + plot_annotation(
  title = paste(species_list[s], ': Theta Over Time'))

ggsave(thing, filename = paste(species_list[s], '_theta.png', sep = ''),
       path = 'C:/Users/Vic/Google Drive/OSU/Thesis/figures/all_FDR_values')

}
