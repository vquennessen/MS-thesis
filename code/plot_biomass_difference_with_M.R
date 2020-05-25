# plot difference in biomass (transient - static) 

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
num_sims <- 2
data_folder <- 'Test'
figures_folder <- 'Test/difference'
###############################################################################

# species to compare
species_list <- c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')
Names <- c('Black Rockfish', 'Cabezon', 'Lingcod', 'Canary Rockfish')

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:6)

# dimensions
TimeT <- Time1 + Time2
CR <- length(Control_rules)
FDR <- length(Final_DRs)
sample_size = num_sims
PD = 0.25
plot_individual_runs = FALSE
Error = 0.05
estimates <- c('True', 'Low', 'High')
ENM = 2

nT <- Time2 + 1
nE <- length(estimates)
nC <- length(Control_rules)
nS <- length(species_list)
nF <- length(Final_DRs)

base <- data.frame(Estimate = rep(estimates, each = nF*nT), 
                   FDR = rep(Final_DRs, times = nE, each = nT), 
                   Time = rep(0:Time2, times = nE*nF),
                   Difference = rep(0, nE*nF*nT), 
                   Lower = rep(0, nE*nF*nT), 
                   Upper = rep(0, nE*nF*nT))

for (s in 1:length(species_list)) {
  
  # load objects
  load(paste('~/Projects/MS-thesis/data/', data_folder, '/', 
             species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
  
  Rec_age <- parameters(species_list[s])[[3]]
  Max_age <- parameters(species_list[s])[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)

  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  ##### relative biomass and median, upper, and lower limits  #####
  
  # pull out sample sims
  B_sample   <- sims_biomass[, , , , indices]
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, MPA*nT*nC*nF*num_sims), 
                       c(MPA, nT, nC, nF, num_sims))
  
  # calculate relative arrays after reserve implementation
  for (cr in 1:nC) {
    for (fdr in 1:nF) {
      for (sim in indices) {
        for (a in 1:MPA) {
          Rel_biomass[a, , cr, fdr, sim] <- B_sample[a, , cr, fdr, sim] / B_sample[a, 1, cr, fdr, sim]
        }
      }
    }
  }
  
  # initialize median, lowerIQR, and upperIQR arrays
  B_medians <- B_lower <- B_upper <- array(rep(NA, MPA*nT*nC*nF), 
                                           c(MPA, nT, nC, nF))
  
  # calculate medians, upper limits, and lower limits  
  for (t in 1:nT) {
    for (cr in 1:nC) {
      for (fdr in 1:nF) {     
        for (a in 1:MPA) {
          B_medians[a, t, cr, fdr] <- median(Rel_biomass[a, t, cr, fdr, ])
          B_lower[a, t, cr, fdr] <- quantile(Rel_biomass[a, t, cr, fdr, ], 0.5 - PD)
          B_upper[a, t, cr, fdr] <- quantile(Rel_biomass[a, t, cr, fdr, ], 0.5 + PD)
        }
      }
    }
  } 
  
  ##### far - initialize and fill in data frame #####
  
  # make copies of base data frame called far, near, and inside
  far <- base; near <- base; inside <- base
  
  for (e in 1:nE) {
    for (fdr in 1:nF) {
      for (t in 1:nT) {
        index <- (e - 1)*nF*nT + (fdr - 1)*nT + t

        far$Difference[index] <- (B_medians[1, t, 2*e, fdr] -
                                    B_medians[1, t, 2*e - 1, fdr]) / B_medians[1, t, 2*e - 1, fdr]
        far$Lower[index] <- (B_lower[1, t, 2*e, fdr] -
                               B_lower[1, t, 2*e - 1, fdr]) / B_lower[1, t, 2*e - 1, fdr]
        far$Upper[index] <- (B_upper[1, t, 2*e, fdr] -
                               B_upper[1, t, 2*e - 1, fdr]) / B_upper[1, t, 2*e - 1, fdr]

        near$Difference[index] <- (B_medians[2, t, 2*e, fdr] -
                                     B_medians[2, t, 2*e - 1, fdr]) / B_medians[2, t, 2*e - 1, fdr]
        near$Lower[index] <- (B_lower[2, t, 2*e, fdr] -
                                B_lower[2, t, 2*e - 1, fdr]) / B_lower[2, t, 2*e - 1, fdr]
        near$Upper[index] <- (B_upper[2, t, 2*e, fdr] -
                                B_upper[2, t, 2*e - 1, fdr]) / B_upper[2, t, 2*e - 1, fdr]

        inside$Difference[index] <- (B_medians[3, t, 2*e, fdr] -
                                       B_medians[3, t, 2*e - 1, fdr]) / B_medians[3, t, 2*e - 1, fdr]
        inside$Lower[index] <- (B_lower[3, t, 2*e, fdr] -
                                  B_lower[3, t, 2*e - 1, fdr]) / B_lower[3, t, 2*e - 1, fdr]
        inside$Upper[index] <- (B_upper[3, t, 2*e, fdr] -
                                  B_upper[3, t, 2*e - 1, fdr]) / B_upper[3, t, 2*e - 1, fdr]

      }
    }
  }
  
  # for (e in 1:nE) {
  #   for (fdr in 1:nF) {
  #     for (t in 1:nT) {
  #       index <- (e - 1)*nF*nT + (fdr - 1)*nT + t
  # 
  #       far$Difference[index] <- (B_medians[1, t, e + 3, fdr] -
  #                                   B_medians[1, t, e, fdr]) / B_medians[1, t, e, fdr]
  #       far$Lower[index] <- (B_lower[1, t, e + 3, fdr] -
  #                              B_lower[1, t, e, fdr]) / B_lower[1, t, e, fdr]
  #       far$Upper[index] <- (B_upper[1, t, e + 3, fdr] -
  #                              B_upper[1, t, e, fdr]) / B_upper[1, t, e, fdr]
  # 
  #       near$Difference[index] <- (B_medians[2, t, e + 3, fdr] -
  #                                    B_medians[2, t, e, fdr]) / B_medians[2, t, e, fdr]
  #       near$Lower[index] <- (B_lower[2, t, e + 3, fdr] -
  #                               B_lower[2, t, e, fdr]) / B_lower[2, t, e, fdr]
  #       near$Upper[index] <- (B_upper[2, t, e + 3, fdr] -
  #                               B_upper[2, t, e, fdr]) / B_upper[2, t, e, fdr]
  # 
  #       inside$Difference[index] <- (B_medians[3, t, e + 3, fdr] -
  #                                      B_medians[3, t, e, fdr]) / B_medians[3, t, e, fdr]
  #       inside$Lower[index] <- (B_lower[3, t, e + 3, fdr] -
  #                                 B_lower[3, t, e, fdr]) / B_lower[3, t, e, fdr]
  #       inside$Upper[index] <- (B_upper[3, t, e + 3, fdr] -
  #                                 B_upper[3, t, e, fdr]) / B_upper[3, t, e, fdr]
  # 
  #     }
  #   }
  # }
  
  ##### plotting parameters #####
  y1 <- -0.15
  y2 <- 2.25
  y1.1 <- y1/10
  y2.1 <- y2/10
  jitter_height1 <- 0.01
  jitter_height2 <- jitter_height1/10
  ##### plot far #####
  
  # FAR <- ggplot(data = far, aes(x = Time, y = Difference,
  #                               color = as.factor(FDR), 
  #                               linetype = as.factor(Estimate))) +
  #   geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  #   geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  #   ggtitle('Far From Reserve') +
  #   ylab('Difference in Relative Biomass') +
  #   theme(legend.position = 'none') +
  #   theme(axis.title.x = element_blank()) +
  #   ylim(y1, y2)
  
  OUTSIDE <- ggplot(data = near, aes(x = Time, y = Difference,
                                color = as.factor(FDR), 
                                linetype = Estimate)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height1)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
    ggtitle('Outside Reserve') +
    ylab('Difference in Relative Biomass') +
    xlab('Years since reserve implementation') +
    theme(legend.position = 'none') +
    ylim(y1, y2)
  
  ##### plot near #####
  
  # NEAR <- ggplot(data = near, aes(x = Time, y = Difference,
  #                                color = as.factor(FDR))) +
  #   geom_line() +
  #   geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
  #   labs(color = 'FDR') +
  #   ggtitle('Near Reserve') +
  #   theme(axis.title.y = element_blank()) +
  #   theme(axis.text.y = element_blank()) +
  #   theme(legend.position = 'none')  +
  #   ylim(y1, y2)
  
  ##### plot inside #####
  
  INSIDE <- ggplot(data = inside, aes(x = Time, y = Difference,
                                      color = as.factor(FDR), 
                                      linetype = Estimate)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height2)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
    labs(color = 'FDR') +
    ggtitle('Inside Reserve') +
    # ylab('Difference in Relative Biomass') +
    theme(axis.title.y = element_blank()) +
    xlab('Years since reserve implementation') +
    # theme(legend.position = c(1.5, 0.4))  +
    ylim(y1.1, y2.1) +
    labs(color = 'FDR', linetype = 'Estimate')
  
  ##### patch all the figures together #####
  # patch <- (FAR + NEAR) / (INSIDE + plot_spacer())
  patch <- OUTSIDE + INSIDE
  thing <- patch + plot_annotation(
    title = paste(Names[s], ': Difference in Biomass', sep = ''))
  
  ggsave(thing, filename = paste('M_biomass_difference_', Names[s], '.png', sep = ''),
         path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                      figures_folder, sep = ''))
  
}
