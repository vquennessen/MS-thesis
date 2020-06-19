# plot total theta for static vs. transient DRs over time

# load any necessary libraries
library(ggplot2)
library(patchwork)
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
num_sims <- 6193
data_folder <- 'Both'
figures_folder <- 'Both/relative'
cluster <- TRUE
png_width <- 7
png_height <- 6
y1 = 0.75
y2 = 1.4
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')

# set variables
A = 5
MPA = 3
Time2 = 20
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules <- 1:6
R0 <- 1e5

# dimensions
PD = 0.25
Error = 0.05
estimates <- c('True', 'Low', 'High')
ENM = 2
types <- c('Static', 'Transient')
sample_size = num_sims

nC <- length(Control_rules)
nT <- Time2 + 1
nE <- length(estimates)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1

theta_df <- data.frame(Species = c(rep(Names[1], each = 2*nE*nF1*nT),
                                   rep(Names[2], each = 2*nE*nF1*nT),
                                   rep(Names[3], each = 2*nE*nF1*nT),
                                   rep(Names[4], each = 2*nE*nF2*nT)), 
                       Estimate = c(rep(estimates, each = 2*nF1*nT), 
                                    rep(estimates, each = 2*nF1*nT), 
                                    rep(estimates, each = 2*nF1*nT), 
                                    rep(estimates, each = 2*nF2*nT)),
                       Type = c(rep(types, each = nF1*nT, times = nE), 
                                rep(types, each = nF1*nT, times = nE), 
                                rep(types, each = nF1*nT, times = nE), 
                                rep(types, each = nF2*nT, times = nE)), 
                       FDR = c(rep(Final_DRs1, each = nT, times = 2*nE), 
                               rep(Final_DRs1, each = nT, times = 2*nE), 
                               rep(Final_DRs1, each = nT, times = 2*nE), 
                               rep(Final_DRs2, each = nT, times = 2*nE)), 
                       Time = c(rep(0:Time2, times = 2*nF1*nE), 
                                rep(0:Time2, times = 2*nF1*nE), 
                                rep(0:Time2, times = 2*nF1*nE), 
                                rep(0:Time2, times = 2*nF2*nE)), 
                       Theta = rep(0, 2*nE*nFall*nT), 
                       Lower = rep(0, 2*nE*nFall*nT), 
                       Upper = rep(0, 2*nE*nFall*nT))

for (s in 1:length(species_list)) {
  
  # load objects
  if (cluster == TRUE) {
    load(paste('~/Documents/MS-thesis/data/', data_folder, '/', 
               species_list[s], '/', num_sims, '_N.Rda', sep = ''))
  } else {
    load(paste('~/Projects/MS-thesis/data/', data_folder, '/', 
               species_list[s], '/', num_sims, '_N.Rda', sep = ''))
  }
  
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
  
  ##### set nF value for species #####
  nF <- ifelse(s == 4, nF2, nF1)
  
  # pull out sample sims and calculate median across num_sims
  A_sample <- array(rep(0, n*nT*nC*nF*num_sims), c(n, nT, nC, nF, num_sims))
  A_medians <- array(rep(0, n*nT*nC*nF), c(n, nT, nC, nF))
  A_lower <- array(rep(0, n*nT*nC*nF), c(n, nT, nC, nF))
  A_upper <- array(rep(0, n*nT*nC*nF), c(n, nT, nC, nF))
  
  for (i in 1:n) {
    A_sample[i, , , , ] <- colSums(sims_N[i, , , , , ])
    for (t in 1:nT) {
      for (cr in 1:nC) {
        for (fdr in 1:nF) {
          A_medians[i, t, cr, fdr] <- median(A_sample[i, t, cr, fdr, ])
          A_lower[a, t, cr, fdr] <- quantile(A_sample[a, t, cr, fdr, ], 0.5 - PD)
          A_upper[a, t, cr, fdr] <- quantile(A_sample[a, t, cr, fdr, ], 0.5 + PD)
        }
      }
    }
  }
  
   for (e in 1:nE) {
     for (type in 1:2) {
      for (fdr in 1:nF) {
        for (t in 1:nT) {
          
          cr <- 2*e - type %% 2
          
          num <- SAD %*% A_medians[, t, cr, fdr]
          num_L <- SAD %*% A_lower[, t, cr, fdr]
          num_U <- SAD %*% A_upper[, t, cr, fdr]
          norm1 <- norm(matrix(SAD), type = 'F')
          norm2 <- norm(matrix(A_medians[, t, cr, fdr]), type = 'F')
          norm2_L <- norm(matrix(A_lower[, t, cr, fdr]), type = 'F')
          norm2_U <- norm(matrix(A_upper[, t, cr, fdr]), type = 'F')
          
          denom <- norm1 * norm2
          denom_L <- norm1 * norm2_L
          denom_U <- norm1 * norm2_U
          
          i <- (s - 1)*2*nE*nF1*nT + (e - 1)*2*nF*nT + (type - 1)*nF*nT + 
            (fdr - 1)*nT + t
          # print(i)
          theta_df$Theta[i] <- as.numeric(acos(num / denom))
          theta_df$Lower[i] <- as.numeric(acos(num_L / denom_L))
          theta_df$Upper[i] <- as.numeric(acos(num_U / denom_U))
          
        }
      }
    }
  }
  
}

dfA <- subset(theta_df, Species == Names[1])
dfB <- subset(theta_df, Species == Names[2])
dfC <- subset(theta_df, Species == Names[3])
dfD <- subset(theta_df, Species == Names[4])

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
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +  
  scale_fill_manual(values = alpha(c(colors), 0.25), guide = FALSE) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  scale_color_manual(values = c("#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle(Names[1]) +
  ylab('Theta') +
  ylim(y1, y2) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

B <- ggplot(data = dfB, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = Type)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +  
  scale_fill_manual(values = alpha(c(colors), 0.25), guide = FALSE) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  scale_color_manual(values = c("#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle(Names[2]) +
  ylim(y1, y2) +
  theme(axis.text = element_blank()) +
  theme(axis.title = element_blank()) +
  theme(legend.position = 'none') 

C <- ggplot(data = dfC, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = Type)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +  
  scale_fill_manual(values = alpha(c(colors), 0.25), guide = FALSE) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  scale_color_manual(values = c("#00BA38", "#00BFC4", "#619CFF", "#F564E3")) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle(Names[3]) +
  xlab('Years since reserve implemented') +
  theme(axis.title.x = element_blank()) +
  ylab('Theta') +
  ylim(y1, y2) +
  theme(legend.position = 'none')

D <- ggplot(data = dfD, aes(x = Time, y = Theta, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = Type)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +  
  scale_fill_manual(values = alpha(c(colors), 0.25), guide = FALSE) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  scale_color_manual(values = c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", 
                                "#619CFF", "#F564E3")) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle(Names[4]) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ylim(y1, y2) +
  xlab('Years since reserve implemented') +
  theme(axis.title.x = element_text(hjust = -4.5)) +
  theme(plot.margin = unit(c(0, 80, 0, 0), 'pt')) +
  labs(linetype = 'Type', color = 'FDR') +
  theme(legend.position = c(1.25, 1.1)) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         size = guide_legend(order = 3))

##### patch all the figures together #####
patch2 <- (A + B) / (C + D)
thing2 <- patch2 + plot_annotation(
  title = 'Distance of total population from stable age distribution')

if (cluster == TRUE) {
  ggsave(thing2, filename = 'M_relative_theta.png',
         path = paste('~/Documents/MS-thesis/figures/', figures_folder, sep = ''),
         width = png_width, height = png_height)
} else {
  ggsave(thing2, filename = 'M_relative_theta.png',
         path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                      figures_folder, sep = ''),
         width = png_width, height = png_height)
}

