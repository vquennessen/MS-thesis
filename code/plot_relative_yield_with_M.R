# plot relative total yield for static and transient DRs with M

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

# set variables
A = 5
MPA = 3
Time2 = 20
Final_DRs <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules <- 1:6

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

yield_df <- data.frame(Species = rep(species_list, each = 2*nE*nF*nT),
                       Type = rep(c('Static', 'Transient'), times = nS,
                                  each = nE*nF*nT),
                       Estimate = rep(estimates, times = nS*2, each = nF*nT),
                       FDR = rep(Final_DRs, times = nS*2*nE, each = nT), 
                       Time = rep(0:Time2, times = nS*2*nE*nF),
                       Value = rep(0, nS*2*nE*nF*nT), 
                       Lower = rep(0, nS*2*nE*nF*nT), 
                       Upper = rep(0, nS*2*nE*nF*nT))

for (s in 1:length(species_list)) {
  
  # load objects
  load(paste('~/Projects/MS-thesis/data/', data_folder, '/', 
             species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
  
  Rec_age <- parameters(species_list[s])[[3]]
  Max_age <- parameters(species_list[s])[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  ##### relative biomass and median, upper, and lower limits  #####
  
  # pull out sample sims
  Y_sample   <- colSums(sims_yield[, , , , indices])
  
  # initialize relative arrays
  Rel_yield <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
  
  # calculate relative arrays after reserve implementation
  for (cr in 1:nC) {
    for (fdr in 1:nF) {
      for (sim in 1:num_sims) {
        Rel_yield[, cr, fdr, sim] <- Y_sample[, cr, fdr, sim] / 
          Y_sample[1, cr, fdr, sim]
      }
    }
  }
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- Y_lower <- Y_upper <- array(rep(NA, nT*nC*nF), c(nT, nC, nF))
  
  # calculate medians, upper limits, and lower limits  
  for (t in 1:nT) {
    for (cr in 1:nC) {
      for (fdr in 1:nF) {        
        Y_medians[t, cr, fdr] <- median(Rel_yield[t, cr, fdr, ])
        Y_lower[t, cr, fdr] <- quantile(Rel_yield[t, cr, fdr, ], 0.5 - PD)
        Y_upper[t, cr, fdr] <- quantile(Rel_yield[t, cr, fdr, ], 0.5 + PD)
      }
    }
  } 
  
  ##### initialize and fill in data frame #####
  
  for (type in 1:2) {
    for (e in 1:nE) {
      for (fdr in 1:nF) {
        for (t in 1:nT) {
          index <- (s - 1)*2*nE*nF*nT + (type - 1)*nE*nF*nT + (e - 1)*nF*nT + 
            (fdr - 1)*nT + t
          yield_df$Value[index] <- Y_medians[t, (type - 1)*3 + e, fdr]
          yield_df$Lower[index] <- Y_lower[t, (type - 1)*3 + e, fdr]
          yield_df$Upper[index] <- Y_upper[t, (type - 1)*3 + e, fdr]
        }
      }
    }
  }
  
}

nI <- nC*nT*nF

dfA <- yield_df[1:nI, ]
dfB <- yield_df[(nI + 1):(2*nI), ]
dfC <- yield_df[(2*nI + 1):(3*nI), ]
dfD <- yield_df[(3*nI + 1):(4*nI), ]

# plotting parameters
y1 = 0.15
y2 = 1.4
jitter_height <- 0.015
size1 <- 1.5
size2 <- 0.75
ltype1 <- 2
ltype2 <- 1

# plot it

A <- ggplot(data = dfA, aes(x = Time, y = Value, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2), guide = 'none') +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') + 
  ggtitle('Black Rockfish') +
  ylab('Relative Yield') +
  ylim(y1, y2) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

B <- ggplot(data = dfB, aes(x = Time, y = Value, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ggtitle('Cabezon') +
  ylim(y1, y2) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = 'none')

C <- ggplot(data = dfC, aes(x = Time, y = Value, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2), guide = 'none') +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') + 
  ggtitle('Lingcod') +
  xlab('Years since reserve implemented') +
  ylab('Relative Yield') +
  ylim(y1, y2) +
  theme(legend.position = 'none')

D <- ggplot(data = dfD, aes(x = Time, y = Value, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_size_manual(values = c(size1, size2)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ggtitle('Copper Rockfish') +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ylim(y1, y2) +
  xlab('Years since reserve implemented') +
  labs(color = 'FDR', linetype = 'Estimate', size = 'Type') +
  theme(legend.position = c(1.2, 1)) +
  theme(plot.margin = unit(c(0, 80, 0, 0), 'pt'))

##### patch all the figures together #####
patch2 <- (A + B) / (C + D)
thing2 <- patch2 + plot_annotation(
  title = 'Relative Yield')

ggsave(thing2, filename = 'M_relative_yield.png',
       path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                    figures_folder, sep = ''))
