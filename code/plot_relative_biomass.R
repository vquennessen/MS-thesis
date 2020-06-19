# plot relative biomass for static vs. transient DRs for each area

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
num_sims <- 6193
data_folder <- 'Both'
figures_folder <- 'Both/relative'
cluster <- TRUE
png_width <- 6
png_height <- 5
y1 <- 0.25
y2 <- 1.5
y1.1 <- 0.95
y2.1 <- 2
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:6)
types <- c('Static', 'Transient')

# dimensions
sample_size = num_sims
PD = 0.25
plot_individual_runs = FALSE
Error = 0.05
estimates <- c('True', 'Low', 'High')

nT <- Time2 + 1
nE <- length(estimates)
nC <- length(Control_rules)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1

base1 <- data.frame(Estimate = rep(estimates, each = 2*nF1*nT),
                    Type = rep(types, each = nF1*nT, times = nE), 
                    FDR = rep(Final_DRs1, each = nT, times = 2*nE), 
                    Time = rep(0:Time2, times = 2*nF1*nE), 
                    Value = rep(0, 2*nE*nF1*nT), 
                    Lower = rep(0, 2*nE*nF1*nT),
                    Upper = rep(0, 2*nE*nF1*nT))

base2 <- data.frame(Estimate = rep(estimates, each = 2*nF2*nT),
                    Type = rep(types, each = nF2*nT, times = nE), 
                    FDR = rep(Final_DRs2, each = nT, times = 2*nE), 
                    Time = rep(0:Time2, times = 2*nF2*nE), 
                    Value = rep(0, 2*nE*nF2*nT), 
                    Lower = rep(0, 2*nE*nF2*nT),
                    Upper = rep(0, 2*nE*nF2*nT))

for (s in 1:length(species_list)) {
    
  # load objects
  if (cluster == TRUE) {
    load(paste('~/Documents/MS-thesis/data/', data_folder, '/', 
               species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
  } else {
  load(paste('~/Projects/MS-thesis/data/', data_folder, '/', 
             species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
  }
  
  Rec_age <- parameters(species_list[s])[[3]]
  Max_age <- parameters(species_list[s])[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  M <- parameters(species_list[s])[[2]]
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  # set nF value for species 
  nF <- ifelse(s == 4, nF2, nF1)  
  
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
  
  ##### fill in data frame #####
  
  # make copies of base data frame called far, near, and inside
  if (s != 4) {far <- base1; near <- base1; inside <- base1
  } else { far <- base2; near <- base2; inside <- base2 }
  
  for (e in 1:nE) { 
    for (type in 1:2) {
      for (fdr in 1:nF) {
        for (t in 1:nT) {
          
          cr <- 2*e - type %% 2
          
          index <- (e - 1)*2*nF*nT + (type - 1)*nF*nT + (fdr - 1)*nT + t
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
  }
  
  ##### plotting parameters #####
  jitter_height <- 0.005
  size1 <- 1.5
  size2 <- 0.75
  
  if (s != 4) {
    colors <- c("#00BA38", "#00BFC4", "#619CFF", "#F564E3")
  } else {
    colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
  }
  
  ##### plot outside (far) #####
  
  OUTSIDE <- ggplot(data = far, aes(x = Time, y = Value, 
                                    color = as.factor(FDR), 
                                    linetype = as.factor(Estimate), 
                                    size = as.factor(Type))) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), colour = NA)) +
    scale_fill_manual(values = alpha(c(colors), 0.05)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_color_manual(values = colors) +
    scale_size_manual(values = c(size1, size2)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    labs(color = 'FDR', linetype = 'Estimate', size = 'Type') +
    ggtitle('Outside Reserve') +
    ylab('Relative Biomass') +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = 'none') +
    ylim(y1, y2)
  
  ##### plot inside #####
  
  INSIDE <- ggplot(data = inside, aes(x = Time, y = Value,
                                      color = as.factor(FDR),
                                      linetype = as.factor(Estimate), 
                                      size = Type)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                    colour = NA), show.legend = FALSE) +
    scale_fill_manual(values = alpha(colors, 0.05)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_color_manual(values = colors) +
    scale_size_manual(values = c(size1, size2)) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    ggtitle('Inside Reserve') +
    theme(axis.title.y = element_blank()) +
    xlab('Years since reserve implementation') +
    theme(axis.title.x = element_text(hjust = 2.5)) +
    ylim(y1.1, y2.1) +
    labs(color = 'FDR', linetype = 'Estimate', size = 'Type') +
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), 
           size = guide_legend(order = 3))
  
  ##### patch all the figures together #####
  # patch <- (FAR + NEAR) / (INSIDE + plot_spacer())
  patch <- OUTSIDE + INSIDE
  thing <- patch + plot_annotation(
    title = paste(Names[s], ': Relative Biomass', sep = ''))
  
  if (cluster == TRUE) {
    
    ggsave(thing, 
           filename = paste( 'M_biomass_relative_', Names[s], '.png', sep = ''),
           path = paste('~/Documents/MS-thesis/figures/', figures_folder, 
                        sep = ''),
           width = png_width, height = png_height)
    
  } else {
    
  ggsave(thing, filename = paste( 'M_biomass_relative_', Names[s], '.png', sep = ''),
         path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                      figures_folder, sep = ''),
         width = png_width, height = png_height)
  }
  
}
