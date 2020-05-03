# plot difference in biomass for static vs. transient DRs for each area

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
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
CR <- length(Control_rules)
FDR <- length(Final_DRs)
sample_size = num_sims
PD = 0.25
plot_individual_runs = FALSE
Error = 0.05
ENM = 2

for (s in 1:length(species_list)) {
  
  # load objects
  load(paste('~/Projects/MS-thesis/data/no_stochasticity/all_FDR_values/', 
             species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
  
  Rec_age <- parameters(species_list[s])[[3]]
  Max_age <- parameters(species_list[s])[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  M <- parameters(species_list[s])[[2]]
  Nat_mortality <- c(M - Error, M, M + Error)
  NM <- length(Nat_mortality)
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  ##### relative biomass and median, upper, and lower limits  #####
  
  # pull out sample sims
  B_sample   <- sims_biomass[, , , , indices]
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, MPA*(Time2 + 1)*CR*FDR*num_sims), 
                       c(MPA, Time2 + 1, CR, FDR, num_sims))
  
  # calculate relative arrays after reserve implementation
  for (cr in 1:CR) {
    for (fdr in 1:FDR) {
      for (sim in indices) {
        for (a in 1:MPA) {
          Rel_biomass[a, , cr, fdr, sim] <- B_sample[a, , cr, fdr, sim] / B_sample[a, 1, cr, fdr, sim]
        }
      }
    }
  }
  
  # initialize median, lowerIQR, and upperIQR arrays
  B_medians <- B_lower <- B_upper <- array(rep(NA, MPA*(Time2 + 1)*CR*FDR), 
                                           c(MPA, Time2 + 1, CR, FDR))
  
  # calculate medians, upper limits, and lower limits  
  for (t in 1:(Time2 + 1)) {
    for (cr in 1:CR) {
      for (fdr in 1:FDR) {     
        for (a in 1:MPA) {
          B_medians[a, t, cr, fdr] <- median(Rel_biomass[a, t, cr, fdr, ])
          B_lower[a, t, cr, fdr] <- quantile(Rel_biomass[a, t, cr, fdr, ], 0.5 - PD)
          B_upper[a, t, cr, fdr] <- quantile(Rel_biomass[a, t, cr, fdr, ], 0.5 + PD)
        }
      }
    }
  } 
  
  ##### plotting parameters #####
  y1 <- -0.1
  y2 <- 0.7
  y1.5 <- -0.01
  y2.5 <- 0.15
  
  ##### far - initialize and fill in data frame #####
  
  far <- data.frame(Estimate = rep(c('Low', 'True', 'High'), 
                                   each = (Time2 + 1)*FDR), 
                    FDR = rep(Final_DRs, each = Time2 + 1, times = NM), 
                    Time = rep(0:Time2, times = NM*FDR),
                    Difference = rep(0, (Time2 + 1)*NM*FDR), 
                    Lower = rep(0, (Time2 + 1)*NM*FDR), 
                    Upper = rep(0, (Time2 + 1)*NM*FDR))
  
  for (nm in 1:NM) {
    for (fdr in 1:FDR) {
      for (t in 1:(Time2 + 1)) {
        index <- (nm - 1)*(Time2 + 1)*FDR + (fdr - 1)*(Time2 + 1) + t
        
        difference <- B_medians[1, t, nm + 3, fdr] - B_medians[1, t, nm, fdr]
        diff_low <- B_lower[1, t, nm + 3, fdr] - B_lower[1, t, nm, fdr]
        diff_up <- B_upper[1, t, nm + 3, fdr] - B_upper[1, t, nm, fdr]
        
        far$Difference[index] <- difference / B_medians[1, t, nm, fdr]
        far$Lower[index] <- diff_low / B_lower[1, t, nm, fdr]
        far$Upper[index] <- diff_up / B_upper[1, t, nm, fdr]
      }
    }
  }
  
  # plot it
  FAR <- ggplot(data = far, aes(x = Time, y = Difference, 
                                color = as.factor(FDR), 
                                linetype = as.factor(Estimate))) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
    labs(color = 'FDR', linetype = 'Estimate of M') +
    ggtitle('Far From Reserve') +
    ylab('Difference in Relative Biomass') +
    theme(legend.position = 'none') +
    theme(axis.title.x = element_blank()) +
    ylim(y1, y2)
  
  ##### near - initialize and fill in data frame #####
  
  near <- data.frame(Estimate = rep(c('Low', 'True', 'High'), 
                                    each = (Time2 + 1)*FDR), 
                     FDR = rep(Final_DRs, each = Time2 + 1, times = NM), 
                     Time = rep(0:Time2, times = NM*FDR),
                     Difference = rep(0, (Time2 + 1)*NM*FDR), 
                     Lower = rep(0, (Time2 + 1)*NM*FDR), 
                     Upper = rep(0, (Time2 + 1)*NM*FDR))
  
  for (nm in 1:NM) {
    for (fdr in 1:FDR) {
      for (t in 1:(Time2 + 1)) {
        index <- (nm - 1)*(Time2 + 1)*FDR + (fdr - 1)*(Time2 + 1) + t
        
        difference <- B_medians[2, t, nm + 3, fdr] - B_medians[2, t, nm, fdr]
        diff_low <- B_lower[2, t, nm + 3, fdr] - B_lower[2, t, nm, fdr]
        diff_up <- B_upper[2, t, nm + 3, fdr] - B_upper[2, t, nm, fdr]

        near$Difference[index] <- difference / B_medians[2, t, nm, fdr]
        near$Lower[index] <- diff_low / B_lower[2, t, nm, fdr]
        near$Upper[index] <- diff_up / B_upper[2, t, nm, fdr]
      }
    }
  }
  
  # plot it
  NEAR <- ggplot(data = near, aes(x = Time, y = Difference, 
                                  color = as.factor(FDR), 
                                  linetype = as.factor(Estimate))) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
    labs(color = 'FDR', linetype = 'Estimate of M') +
    ggtitle('Near Reserve') +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(legend.position = 'none')  +
    ylim(y1, y2)
  
  ##### inside - initialize and fill in data frame #####
  
  inside <- data.frame(Estimate = rep(c('Low', 'True', 'High'), 
                                      each = (Time2 + 1)*FDR), 
                       FDR = rep(Final_DRs, each = Time2 + 1, times = NM), 
                       Time = rep(0:Time2, times = NM*FDR),
                       Difference = rep(0, (Time2 + 1)*NM*FDR), 
                       Lower = rep(0, (Time2 + 1)*NM*FDR), 
                       Upper = rep(0, (Time2 + 1)*NM*FDR))
  
  for (nm in 1:NM) {
    for (fdr in 1:FDR) {
      for (t in 1:(Time2 + 1)) {
        index <- (nm - 1)*(Time2 + 1)*FDR + (fdr - 1)*(Time2 + 1) + t
        difference <- B_medians[3, t, nm + 3, fdr] - B_medians[3, t, nm, fdr]
        diff_low <- B_lower[3, t, nm + 3, fdr] - B_lower[3, t, nm, fdr]
        diff_up <- B_upper[3, t, nm + 3, fdr] - B_upper[3, t, nm, fdr]
        
        inside$Difference[index] <- difference / B_medians[3, t, nm, fdr]        
        inside$Lower[index] <- diff_low / B_lower[3, t, nm, fdr]
        inside$Upper[index] <- diff_up / B_upper[3, t, nm, fdr]
      }
    }
  }
  
  # plot it
  INSIDE <- ggplot(data = inside, aes(x = Time, y = Difference, 
                                      color = as.factor(FDR), 
                                      linetype = as.factor(Estimate))) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
    ggtitle('Inside Reserve') +
    ylab('Difference in Relative Biomass') +
    xlab('Years since reserve implementation') +
    theme(legend.position = c(1.5, 0.4), legend.box = 'horizontal')  +
    ylim(y1.5, y2.5) +
    labs(color = 'FDR', linetype = 'M')
  
  ##### patch all the figures together #####
  patch <- (FAR + NEAR) / (INSIDE + plot_spacer())
  thing <- patch + plot_annotation(
    title = paste(species_list[s], 
                  ': Relative Biomass (Transient - Static DRCR)'))
  
  ggsave(thing, filename = paste(species_list[s], '_biomass.png', sep = ''),
         path = 'C:/Users/Vic/Google Drive/OSU/Thesis/figures/all_FDR_values')
  
}
