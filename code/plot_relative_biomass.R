# plot relative biomass for static vs. transient DRs for each area

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
data_folder <- 'no_stochasticity/all_FDR_values'
folder <- 'all_FDR_values/biomass/without_M'
Names <- c('Black Rockfish', 'Cabezon', 'Lingcod', 'Canary Rockfish')
###############################################################################

# species to compare
species_list <- c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')
num_sims <- 2

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:2)

# dimensions
TimeT <- Time1 + Time2
CR <- length(Control_rules)
FDR <- length(Final_DRs)
sample_size = num_sims
PD = 0.25
plot_individual_runs = FALSE

nS <- length(species_list)
nC <- 2
nF <- FDR
nT <- Time2 + 1

base <- data.frame(Type = rep(c('Static', 'Transient'), each = nF*nT),
                   FDR = rep(Final_DRs, each = nT, times = nC), 
                   Time = rep(0:Time2, times = nC*nF),
                   Value = rep(0, nC*nT*nF), 
                   Lower = rep(0, nC*nT*nF), 
                   Upper = rep(0, nC*nT*nF))

for (s in 1:length(species_list)) {
  
  # load objects
  load(paste('~/Projects/MS-thesis/data/', data_folder, '/', 
             species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
  
  Rec_age <- parameters(species_list[s])[[3]]
  Max_age <- parameters(species_list[s])[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  M <- parameters(species_list[s])[[2]]

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
  
  for (cr in 1:nC) {
    for (fdr in 1:nF) {
      for (t in 1:nT) {
        index <- (cr - 1)*nF*nT + (fdr - 1)*nT + t

        far$Value[index] <- B_medians[1, t, cr, fdr]
        far$Lower[index] <- B_lower[1, t, cr, fdr]
        far$Upper[index] <- B_upper[1, t, cr, fdr]
        near$Value[index] <- B_medians[2, t, cr, fdr]
        near$Lower[index] <- B_lower[2, t, cr, fdr]
        near$Upper[index] <- B_upper[2, t, cr, fdr]
        inside$Value[index] <- B_medians[3, t, cr, fdr]
        inside$Lower[index] <- B_lower[3, t, cr, fdr]
        inside$Upper[index] <- B_upper[3, t, cr, fdr]
        
      }
    }
  }
  
  ##### plotting parameters #####
  y1 <- 0.35
  y2 <- 1.8
  jitter_height <- 0.01
  ltype1 <- 2
  ltype2 <- 1
  
  ##### plot far #####
  
  FAR <- ggplot(data = far, aes(x = Time, y = Value,
                                color = as.factor(FDR),
                                linetype = as.factor(Type))) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_linetype_manual(values = c(ltype1, ltype2)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    labs(color = 'FDR', linetype = 'Type') +
    ggtitle('Far From Reserve') +
    ylab('Relative Biomass') +
    theme(legend.position = 'none') +
    theme(axis.title.x = element_blank()) +
    ylim(y1, y2)
  
  OUTSIDE <- ggplot(data = far, aes(x = Time, y = Value,
                                    color = as.factor(FDR),
                                    linetype = as.factor(Type))) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_linetype_manual(values = c(ltype1, ltype2)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    labs(color = 'FDR', linetype = 'Type') +
    ggtitle('Outside Reserve') +
    ylab('Relative Biomass') +
    xlab('Years since reserve implementation') +
    theme(legend.position = 'none') +
    ylim(y1, y2)
  
  ##### plot near #####
  
  NEAR <- ggplot(data = near, aes(x = Time, y = Value,
                                  color = as.factor(FDR),
                                  linetype = as.factor(Type))) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_linetype_manual(values = c(ltype1, ltype2)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    labs(color = 'FDR', linetype = 'Type') +
    ggtitle('Near Reserve') +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(legend.position = 'none')  +
    ylim(y1, y2)
  
  ##### plot inside #####
  
  INSIDE <- ggplot(data = inside, aes(x = Time, y = Value,
                                      color = as.factor(FDR),
                                      linetype = as.factor(Type))) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_linetype_manual(values = c(ltype1, ltype2)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    ggtitle('Inside Reserve') +
    theme(axis.title.y = element_blank()) +
    xlab('Years since reserve implementation') +
    ylim(y1, y2) +
    labs(color = 'FDR', linetype = 'Type')
  
  ##### patch all the figures together #####
  # patch <- (FAR + NEAR) / (INSIDE + plot_spacer())
  patch <- OUTSIDE + INSIDE
  thing <- patch + plot_annotation(
    title = paste(Names[s], ': Relative Biomass', sep = ''))
  
  ggsave(thing, filename = paste(species_list[s], '_biomass.png', sep = ''),
         path = paste('C:/Users/Vic/Google Drive/OSU/Thesis/figures/', folder, 
                      sep = ''))
  
}
