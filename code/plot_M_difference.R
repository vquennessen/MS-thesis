# plot differences in biomass, yield, and effort by species

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(viridisLite)

###############################################################################
# CHECK THESE EVERY TIME
folder <- 'None'
cluster <- FALSE
Years <- c(10, 20)
png_width <- 4
png_height <- 4
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')

# determine num_sims based on data folder
num_sims <- ifelse(folder == 'None', 3, 
                   ifelse(folder == 'Both', 6193, 5000))

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
metrics <- c('Biomass', 'Yield')

nY <- length(Years)
nE <- length(estimates)
nC <- length(Control_rules)
nS <- length(species_list)
s1 <- 3; s2 <- 1
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1

DF1 <- data.frame(Metric = rep(metrics, each = nE*nF1),
                  Estimate = rep(estimates, each = nF1, times = 2),
                  FDR = rep(Final_DRs1, times = 2*nE), 
                  Value = rep(NA, 2*nE*nF1))

DF2 <- data.frame(Metric = rep(metrics, each = nE*nF2),
                  Estimate = rep(estimates, each = nF2, times = 2),
                  FDR = rep(Final_DRs2, times = 2*nE), 
                  Value = rep(NA, 2*nE*nF2))

for (y in 1:nY) {
  
  for (s in 1:length(Names)) {
    
    # load biomass and yield files
    if (cluster == TRUE) {
      load(paste('~/Documents/MS-thesis/data/', folder, '/', 
                 species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
      load(paste('~/Documents/MS-thesis/data/', folder, '/', 
                 species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
      
    } else {
      load(paste('~/Projects/MS-thesis/data/', folder, '/', 
                 species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
      load(paste('~/Projects/MS-thesis/data/', folder, '/', 
                 species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
    }
    
    # set nF value for species 
    nF <- ifelse(s == 4, nF2, nF1)  
    
    ##### relative biomass and median, upper, and lower limits  #####
    
    # pull out sample sims as sums across all areas for particular years
    B_sample   <- colSums(sims_biomass) 
    Y_sample <- sims_yield

    # initialize relative arrays
    Rel_biomass <- array(rep(0, nC*nF*num_sims), c(nC, nF, num_sims))
    Rel_yield <- array(rep(0, nC*nF*num_sims), c(nC, nF, num_sims))
    
    # calculate relative arrays after reserve implementation
    for (cr in 1:nC) {
      for (fdr in 1:nF) {
        for (sim in 1:num_sims) {
          Rel_biomass[cr, fdr, sim] <- B_sample[Years[y] + 1, cr, fdr, sim] / 
            B_sample[1, cr, fdr, sim]
          Rel_yield[cr, fdr, sim] <- Y_sample[Years[y] + 1, cr, fdr, sim] / 
            Y_sample[1, cr, fdr, sim]
        }
      }
    }
    
    ##### initialize difference array #####
    Difference_B <- array(rep(0, nE*nF*num_sims), c(nE, nF, num_sims))
    Difference_Y <- array(rep(0, nE*nF*num_sims), c(nE, nF, num_sims))
    
    ##### calculate differences between transient and static DRs #####
    for (e in 1:nE) {
      for (fdr in 1:nF) {     
        for (sim in 1:num_sims) {
          # biomass
          diff_B <- Rel_biomass[2*e, fdr, sim] - Rel_biomass[2*e - 1, fdr, sim]
          Difference_B[e, fdr, sim] <- diff_B / Rel_biomass[2*e - 1, fdr, sim]
          # yield
          diff_Y <- Rel_yield[2*e, fdr, sim] - Rel_yield[2*e - 1, fdr, sim]
          Difference_Y[e, fdr, sim] <- diff_Y / Rel_yield[2*e - 1, fdr, sim]
        }
      }
    }
    
    # DF based on species
    if (s == 4) { DF <- DF2 } else { DF <- DF1 }
    
    ##### fill in data frames with median and quantile values #####
    for (e in 1:nE) {
      for (fdr in 1:nF) {
        index <- (e - 1)*nF + fdr 
        all_biomass <- nrow(subset(DF, Metric == 'Biomass'))
        # relative biomass
        DF$Value[index] <- median(Difference_B[e, fdr, ])
        # relative yield
        DF$Value[index + all_biomass] <- median(Difference_Y[e, fdr, ])
      }
    }
    
    # make FDR and Estimate factor variables
    DF$FDR <- factor(DF$FDR)
    DF$Estimate <- factor(DF$Estimate, levels = c('Low', 'True', 'High'))
    
    # plotting parameters
    if (s == 4) {
      colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", 
                  "#619CFF", "#F564E3")
    } else { colors <- c("#00BA38", "#00BFC4", "#619CFF", "#F564E3") }
    
    # plot difference results
    thing1 <- ggplot(DF, aes(x = Estimate, y = Value, color = FDR, 
                             shape = Metric)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_point(position = position_jitter(w = 0.3, h = 0), size = 3) +
      scale_shape_manual(values = c(1, 3)) +
      scale_color_manual(values = colors) +
      ylab('Median Difference') +
      xlab('Estimate of Natural Mortality (M)') +
      ggtitle(Names[s]) + 
      guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))
    
    # save results to figures folder
    if (cluster == TRUE) {
      ggsave(thing1, filename = paste(Names[s], '_year', Years[y], 
                                      '_difference.png', sep = ''),
             path = paste('~/Documents/MS-thesis/figures/', folder, 
                          '/difference', sep = ''),
             width = png_width, height = png_height)
    
      } else {
      ggsave(thing1, filename = paste(Names[s], '_year', Years[y], 
                                      '_difference.png', sep = ''),
             path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                          folder, '/difference', sep = ''),
             width = png_width, height = png_height)
    }
    
  }
  
}