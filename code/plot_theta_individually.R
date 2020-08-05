# plot total theta for static vs. transient DRs over time

# load any necessary libraries
library(ggplot2)
library(patchwork)
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
folder <- 'None'
cluster <- FALSE
png_width <- 5
png_height <- 4
with_M <- FALSE
FDRs1 <- c(0.6, 0.9)
FDRs2 <- c(0.4, 0.9)
###############################################################################

# determine num_sims based on data folder
num_sims <- ifelse(folder == 'None', 3, 
                   ifelse(folder == 'Both', 6193, 5000))

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

nC <- length(Control_rules)
nT <- Time2 + 1
nE <- length(estimates)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1

base1 <- data.frame(Type = rep(types, each = nE*nF1*nT),
                    Estimate = rep(estimates, each = nF1*nT, times = 2),
                    FDR = rep(Final_DRs1, each = nT, times = 2*nE),
                    Time = rep(0:Time2, times = 2*nF1*nE),
                    Theta = rep(0, 2*nE*nF1*nT), 
                    Lower = rep(0, 2*nE*nF1*nT),
                    Upper = rep(0, 2*nE*nF1*nT))

base2 <- data.frame(Type = rep(types, each = nE*nF2*nT),
                    Estimate = rep(estimates, each = nF2*nT, times = 2),
                    FDR = rep(Final_DRs2, each = nT, times = 2*nE),
                    Time = rep(0:Time2, times = 2*nF2*nE),
                    Theta = rep(0, 2*nE*nF2*nT), 
                    Lower = rep(0, 2*nE*nF2*nT),
                    Upper = rep(0, 2*nE*nF2*nT))

for (s in 1:length(species_list)) {
  
  # load objects
  if (cluster == TRUE) {
    load(paste('~/Documents/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_N.Rda', sep = ''))
  } else {
    load(paste('~/Projects/MS-thesis/data/', folder, '/', 
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
    A_sample[i, , , , ] <- sims_N[i, MPA, , , , ]
    for (t in 1:nT) {
      for (cr in 1:nC) {
        for (fdr in 1:nF) {
          A_medians[i, t, cr, fdr] <- median(A_sample[i, t, cr, fdr, ])
          A_lower[i, t, cr, fdr] <- quantile(A_sample[i, t, cr, fdr, ], 0.5 - PD)
          A_upper[i, t, cr, fdr] <- quantile(A_sample[i, t, cr, fdr, ], 0.5 + PD)
        }
      }
    }
  }
  
  # set which base data frame to use based on species
  if (s == 4) {theta_df <- base2} else {theta_df <- base1}
  
  for (ty in 1:2) {
    for (e in 1:nE) {
      for (fdr in 1:nF) {
        for (t in 1:nT) {
          
          cr <- 2*e - ty %% 2
          
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
          
          i <- (ty - 1)*nE*nF*nT + (e - 1)*nF*nT + (fdr - 1)*nT + t
          # print(i)
          theta_df$Theta[i] <- as.numeric(acos(num / denom))
          theta_df$Lower[i] <- as.numeric(acos(num_L / denom_L))
          theta_df$Upper[i] <- as.numeric(acos(num_U / denom_U))
          
        }
      }
    }
  }
  
  # re-order estimates of natural mortality so they're not forced into 
  # alphabetical order
  theta_df$Estimate <- factor(theta_df$Estimate, 
                              levels = c('Low', 'True', 'High'))
  
  # plotting parameters
  jitter_height <- 0.001
  size1 <- 1.25
  size2 <- 0.5
  
  if (s == 4) {
    colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
  } else { colors <- c("#00BA38", "#00BFC4", "#619CFF", "#F564E3") }
  
  # # set FDRs for minimum and maximum biomass and yields and extract indices to 
  # # get colors
  # if (s == 4) { 
  #   FDRs <- FDRs2
  #   ind <- which(Final_DRs2 %in% FDRs)
  # } else { 
  #   FDRs <- FDRs1
  #   ind <- which(Final_DRs1 %in% FDRs)
  # }
  # 
  # # pull colors out for scenarios that are not 'None'
  # if (folder != 'None') {new_colors <- colors[ind]} else {new_colors <- colors}
  
  # take out M or keep it
  if (with_M == FALSE) { 
    
    theta_df <- subset(theta_df, Estimate == 'True')
    
    # plot it with only true estimate of M
    thing1 <- ggplot(data = theta_df, 
                     aes(x = Time, y = Theta, color = as.factor(FDR), 
                         linetype = as.factor(Type))) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                      colour = NA), show.legend = FALSE) +  
      scale_fill_manual(values = alpha(c(colors), 0.25), guide = FALSE) +
      geom_line(position = position_jitter(w = 0, h = jitter_height)) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
      ylab('Distance from stable age distribution') +
      xlab('Years since reserve implemented') +
      theme(plot.margin = unit(c(0, 80, 0, 0), 'pt')) +
      labs(linetype = 'Type', color = expression('D'[final])) +
      theme(legend.position = c(1.15, 0.5)) +
      guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), 
             size = guide_legend(order = 3)) 
    
    if (cluster == TRUE) {
      ggsave(thing1, filename = paste(Names[s], '_theta.png', sep = ''),
             path = paste('~/Documents/MS-thesis/figures/', folder, sep = ''),
             width = png_width, height = png_height)
    } else {
      ggsave(thing1, filename = paste(Names[s], '_theta.png', sep = ''),
             path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                          folder, sep = ''),
             width = png_width, height = png_height)
    }
    
  } else {
    
    # plot it with estimates of M
    thing1 <- ggplot(data = theta_df, aes(x = Time, y = Theta, color = as.factor(FDR), 
                                          linetype = as.factor(Estimate), size = Type)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                      colour = NA), show.legend = FALSE) +  
      scale_fill_manual(values = alpha(c(colors), 0.25), guide = FALSE) +
      geom_line(position = position_jitter(w = 0, h = jitter_height)) +
      scale_size_manual(values = c(size1, size2)) +
      scale_color_manual(values = colors) +
      geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
      ylab('Distance from stable age distribution') +
      xlab('Years since reserve implemented') +
      theme(plot.margin = unit(c(20, 80, 0, 0), 'pt')) +
      labs(linetype = 'Type', color = expression('D'[final])) +
      theme(legend.position = c(1.15, 0.5)) +
      guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), 
             size = guide_legend(order = 3))
    
    if (cluster == TRUE) {
      ggsave(thing1, filename = paste(Names[s], '_M_theta.png', sep = ''),
             path = paste('~/Documents/MS-thesis/figures/', folder, sep = ''),
             width = png_width, height = png_height)
    } else {
      ggsave(thing1, filename = paste(Names[s], '_M_theta.png', sep = ''),
             path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                          folder, sep = ''),
             width = png_width, height = png_height)
    }
    
  }
  
  
  
}