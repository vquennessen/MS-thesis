# plot difference in biomass (transient - static) 

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
figures_folder <- 'Both/difference'
cluster <- TRUE
png_width <- 5
png_height <- 4
y1 <- -1
y2 <- 1.5
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

# dimensions
sample_size = num_sims
PD = 0.25
Error = 0.05
estimates <- c('True', 'Low', 'High')

nT <- Time2 + 1
nE <- length(estimates)
nC <- length(Control_rules)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1

base1 <- data.frame(Estimate = rep(estimates, each = nF1*nT), 
                    FDR = rep(Final_DRs1, each = nT, times = nE), 
                    Time = rep(0:Time2, times = nF1*nE), 
                    Difference = rep(0, nE*nF1*nT), 
                    Lower = rep(0, nE*nF1*nT), 
                    Upper = rep(0, nE*nF1*nT))

base2 <- data.frame(Estimate = rep(estimates, each = nF2*nT), 
                    FDR = rep(Final_DRs2, each = nT, times = nE), 
                    Time = rep(0:Time2, times = nF2*nE), 
                    Difference = rep(0, nE*nF2*nT), 
                    Lower = rep(0, nE*nF2*nT), 
                    Upper = rep(0, nE*nF2*nT))

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
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  # set nF value for species 
  nF <- ifelse(s == 4, nF2, nF1)  
  
  ##### relative biomass and median, upper, and lower limits  #####
  
  # pull out sample sims
  B_sample   <- sims_biomass[, , , , indices]
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, MPA*nT*nC*nF*num_sims), c(MPA, nT, nC, nF, num_sims))
  
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
  
  ##### initialize difference array #####
  Difference <- array(rep(0, MPA*nT*nE*nF*num_sims), c(MPA, nT, nE, nF, num_sims))
  
  ##### calculate differences between transient and static DRs #####
  for (a in 1:MPA) {
    for (t in 1:nT) {
      for (e in 1:nE) {
        for (fdr in 1:nF) {     
          for (sim in 1:num_sims) {
            diff <- Rel_biomass[a, t, 2*e, fdr, sim] - Rel_biomass[a, t, 2*e - 1, fdr, sim]
            Difference[a, t, e, fdr, sim] <- diff / Rel_biomass[a, t, 2*e - 1, fdr, sim]
          }
        }
      }
    }
  }
  
  # make copies of base data frame called far, near, and inside
  if (s != 4) {far <- base1; near <- base1; inside <- base1
  } else { far <- base2; near <- base2; inside <- base2 }
  
  ##### fill in data frames with median and quantile values #####
  for (t in 1:nT) {
    for (e in 1:nE) {
      for (fdr in 1:nF) {
        index <- (e - 1)*nF*nT + (fdr - 1)*nT + t
        
        far$Difference[index] <- median(Difference[1, t, e, fdr, ])
        far$Lower[index] <- quantile(Difference[1, t, e, fdr, ], 0.5 - PD)
        far$Upper[index] <- quantile(Difference[1, t, e, fdr, ], 0.5 + PD)
        
        near$Difference[index] <- median(Difference[2, t, e, fdr, ])
        near$Lower[index] <- quantile(Difference[2, t, e, fdr, ], 0.5 - PD)
        near$Upper[index] <- quantile(Difference[2, t, e, fdr, ], 0.5 + PD)
        
        inside$Difference[index] <- median(Difference[3, t, e, fdr, ])
        inside$Lower[index] <- quantile(Difference[3, t, e, fdr, ], 0.5 - PD)
        inside$Upper[index] <- quantile(Difference[3, t, e, fdr, ], 0.5 + PD)
      }
    }
  }
  
  ##### plotting parameters #####
  y1.1 <- y1/10
  y2.1 <- y2/10
  jitter_height1 <- 0.01
  jitter_height2 <- jitter_height1/10
  if (s != 4) {
    colors <- c("#00BA38", "#00BFC4", "#619CFF", "#F564E3")
  } else {
    colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
  }
  
  ##### plot outside (far) #####
  # to plot near instead of far, set data = near
  
  OUTSIDE <- ggplot(data = far, aes(x = Time, y = Difference,
                                     color = as.factor(FDR), 
                                     linetype = Estimate)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                    colour = NA), show.legend = FALSE) +
    scale_fill_manual(values = alpha(c(colors), 0.15)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height2)) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
    ggtitle('Outside Reserve') +
    ylab('Difference in Relative Biomass') +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = 'none') +
    ylim(y1, y2)
  
  ##### plot inside #####
  
  INSIDE <- ggplot(data = inside, aes(x = Time, y = Difference,
                                      color = as.factor(FDR), 
                                      linetype = Estimate)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                    colour = NA), show.legend = FALSE) +
    scale_fill_manual(values = alpha(c(colors), 0.15)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height2)) +
    scale_color_manual(values = colors) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') +
    labs(color = 'FDR') +
    ggtitle('Inside Reserve') +
    theme(axis.title.y = element_blank()) +
    xlab('Years since reserve implementation') +
    theme(plot.margin = unit(c(0, 60, 0, 0), 'pt')) +
    theme(legend.position = c(1.3, 0.5))  +
    theme(axis.title.x = element_text(hjust = 1.65)) +
    ylim(y1.1, y2.1) +
    guides(fill = FALSE) + 
    labs(color = 'FDR', linetype = 'Estimate') +
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))
  
  ##### patch all the figures together #####
  # patch <- (FAR + NEAR) / (INSIDE + plot_spacer())
  patch <- OUTSIDE + INSIDE
  thing <- patch + plot_annotation(
    title = paste(Names[s], ': Difference in Biomass', sep = ''))
  
  if (cluster == TRUE) {
    
    ggsave(thing, 
           filename = paste( 'M_biomass_difference_', Names[s], '.png', sep = ''),
           path = paste('~/Documents/MS-thesis/figures/', figures_folder, 
                        sep = ''),
           width = png_width, height = png_height)
    
  } else {
    
    ggsave(thing, filename = paste( 'M_biomass_difference_', Names[s], '.png', sep = ''),
           path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                        figures_folder, sep = ''),
           width = png_width, height = png_height)
  }
  
}
