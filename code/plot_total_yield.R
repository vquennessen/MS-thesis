# plot difference in total yield for static vs. transient DRs

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
CR <- length(Control_rules)
FDR <- length(Final_DRs)
sample_size = num_sims
PD = 0.25
plot_individual_runs = FALSE
Error = 0.05
estimates <- c('Low', 'True', 'High')
ENM = 2

nT <- Time2 + 1
nE <- length(estimates)
nS <- length(species_list)
nF <- FDR

yield_df <- data.frame(Species = rep(species_list, each = nT*nF*nE),
                       Estimate = rep(estimates, each = nT*nF, times = nS), 
                       FDR = rep(Final_DRs, each = nT, times = nE*nS), 
                       Time = rep(0:Time2, times = nE*nF*nS),
                       Difference = rep(0, nT*nE*nF*nS), 
                       Lower = rep(0, nT*nE*nF*nS), 
                       Upper = rep(0, nT*nE*nF*nS))

for (s in 1:length(species_list)) {
  
  # load objects
  load(paste('~/Projects/MS-thesis/data/no_stochasticity/all_FDR_values/', 
             species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
  
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
  Y_sample   <- colSums(sims_yield[, , , , indices])
  
  # initialize relative arrays
  Rel_yield <- array(rep(0, (Time2 + 1)*CR*FDR*num_sims), 
                     c(Time2 + 1, CR, FDR, num_sims))
  
  # calculate relative arrays after reserve implementation
  for (cr in 1:CR) {
    for (fdr in 1:FDR) {
      for (sim in indices) {
        Rel_yield[, cr, fdr, sim] <- Y_sample[, cr, fdr, sim] / 
          Y_sample[1, cr, fdr, sim]
      }
    }
  }
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- Y_lower <- Y_upper <- array(rep(NA, (Time2 + 1)*CR*FDR), 
                                           c(Time2 + 1, CR, FDR))
  
  # calculate medians, upper limits, and lower limits  
  for (t in 1:(Time2 + 1)) {
    for (cr in 1:CR) {
      for (fdr in 1:FDR) {        
        Y_medians[t, cr, fdr] <- median(Rel_yield[t, cr, fdr, ])
        Y_lower[t, cr, fdr] <- quantile(Rel_yield[t, cr, fdr, ], 0.5 - PD)
        Y_upper[t, cr, fdr] <- quantile(Rel_yield[t, cr, fdr, ], 0.5 + PD)
      }
    }
  } 
  
  ##### initialize and fill in data frame #####
  
  for (nm in 1:NM) {
    for (fdr in 1:FDR) {
      for (t in 1:(Time2 + 1)) {
        index <- (s-1)*nE*nT*nF + (nm - 1)*nT*nF + (fdr - 1)*nT + t
        difference <- Y_medians[t, nm + 3, fdr] - Y_medians[t, nm, fdr]
        yield_df$Difference[index] <- difference / Y_medians[t, nm, fdr]
        yield_df$Lower[index] <- Y_lower[t, nm + 3, fdr] - Y_lower[t, nm, fdr]
        yield_df$Upper[index] <- Y_upper[t, nm + 3, fdr] - Y_upper[t, nm, fdr]
      }
    }
  }

}

nI <- nT*nE*nF

dfA <- yield_df[1:nI, ]
dfB <- yield_df[(nI + 1):(2*nI), ]
dfC <- yield_df[(2*nI + 1):(3*nI), ]
dfD <- yield_df[(3*nI + 1):(4*nI), ]

# plotting parameters
y1 = -0.5
y2 = 1.5
jitter_height <- 0.02

# plot it
A <- ggplot(data = dfA, aes(x = Time, y = Difference, 
                           color = as.factor(FDR), 
                           linetype = as.factor(Estimate))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  labs(color = 'FDR', linetype = 'Estimate of M') +
  ggtitle('Black Rockfish') +
  ylab('Difference') +
  ylim(y1, y2) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

B <- ggplot(data = dfB, aes(x = Time, y = Difference, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  labs(color = 'FDR', linetype = 'Estimate of M') +
  ggtitle('Cabezon') +
  ylim(y1, y2) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = 'none')

C <- ggplot(data = dfC, aes(x = Time, y = Difference, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  labs(color = 'FDR', linetype = 'Estimate of M') +
  ggtitle('Lingcod') +
  xlab('Years since reserve implemented') +
  ylab('Difference') +
  ylim(y1, y2) +
  theme(legend.position = 'none')

D <- ggplot(data = dfD, aes(x = Time, y = Difference, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle('Copper Rockfish') +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ylim(y1, y2) +
  xlab('Years since reserve implemented') +
  theme(legend.position = c(1.15, 1)) +
  theme(plot.margin = unit(c(0, 80, 0, 0), 'pt')) +
  labs(color = 'FDR', linetype = 'Estimate \n of M') 

##### patch all the figures together #####
patch2 <- (A + B) / (C + D)
thing2 <- patch2 + plot_annotation(
  title = 'Relative Yield (Transient - Static DRCR)')

ggsave(thing2, filename = 'all_species_yield.png',
       path = 'C:/Users/Vic/Google Drive/OSU/Thesis/figures/all_FDR_values')
