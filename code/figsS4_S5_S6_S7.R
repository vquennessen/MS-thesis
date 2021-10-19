# plot M panels - both relative and difference - in single multi-panel figure

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
years_to_plot <- c(10, 20)
png_width <- 8
png_height <- 6
###############################################################################

# species to compare
# species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
# Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')

species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')
titles <- c('figS4_', 'figS5_', 'figS6_', 'figS7_')

# determine num_sims based on data folder
num_sims <- 1

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:6)
Years <- 1:20

# dimensions
sample_size = num_sims
PD = 0.25
Error = 0.05
types <- c('Static', 'Transient')
estimates <- c('True', 'Low', 'High')
metrics <- c('Biomass', 'Yield')

nTy <- length(types)
nY <- length(Years)
nM <- length(metrics)
nE <- length(estimates)
nC <- length(Control_rules)
nS <- length(species_list)
s1 <- 3; s2 <- 1
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)

# pre-allocate df for actual cumulative values
DF1 <- data.frame(Metric = rep(metrics, each = nTy*nE*nF1*nY),
                      Type = rep(types, times = nM, each = nE*nF1*nY),
                      Estimate = rep(estimates, times = nM*nTy, each = nF1*nY), 
                      FDR = rep(Final_DRs1, times = nM*nTy*nE, each = nY),
                      Year = rep(Years, times = nM*nTy*nE*nF1),
                      Value = rep(NA, nM*nTy*nE*nF1*nY))

DF2 <- data.frame(Metric = rep(metrics, each = nTy*nE*nF2*nY),
                      Type = rep(types, times = nM, each = nE*nF2*nY),
                      Estimate = rep(estimates, times = nM*nTy, each = nF2*nY), 
                      FDR = rep(Final_DRs2, times = nM*nTy*nE, each = nY),
                      Year = rep(Years, times = nM*nTy*nE*nF2),
                      Value = rep(NA, nM*nTy*nE*nF2*nY))

# pre-allocate df for difference in cumulative values
diff1 <- data.frame(Metric = rep(metrics, each = nE*nF1*nY),
                    Estimate = rep(estimates, times = nM, each = nF1*nY), 
                    FDR = rep(Final_DRs1, times = nM*nE, each = nY),
                    Year = rep(Years, times = nM*nE*nF1),
                    Value = rep(NA, nM*nE*nF1*nY))

diff2 <- data.frame(Metric = rep(metrics, each = nE*nF2*nY),
                    Estimate = rep(estimates, times = nM, each = nF2*nY), 
                    FDR = rep(Final_DRs2, times = nM*nE, each = nY),
                    Year = rep(Years, times = nM*nE*nF2),
                    Value = rep(NA, nM*nE*nF2*nY))

for (s in 1:length(Names)) {
  
  # load biomass and yield files
  load(paste('~/Projects/MS-thesis/data/', folder, '/', 
             species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
  load(paste('~/Projects/MS-thesis/data/', folder, '/', 
             species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
  
  # set nF value for species 
  nF <- ifelse(s == 4, nF2, nF1)  
  
  ##### relative biomass and median, upper, and lower limits  #####
  
  # pull out sample sims as sums across all areas for particular years
  B_sample <- colSums(sims_biomass) 
  Y_sample <- sims_yield
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, nY*nC*nF), c(nY, nC, nF))
  Rel_yield <- array(rep(0, nY*nC*nF), c(nY, nC, nF))
  
  # calculate relative arrays after reserve implementation
  for (y in 2:nY) {
    for (cr in 1:nC) {
      for (fdr in 1:nF) {
        Rel_biomass[y, cr, fdr] <- B_sample[y, cr, fdr, 1] / 
          B_sample[1, cr, fdr, 1]
        Rel_yield[y, cr, fdr] <- Y_sample[y, cr, fdr, 1] / 
          Y_sample[1, cr, fdr, 1]
      }
    }
  }
  
  # actual DF based on species
  if (s == 4) { DF <- DF2 } else { DF <- DF1 }
  
  ##### fill in data frames with cumulative values #####
  for (ty in 1:nTy) {
    for (e in 1:nE) {
      
      # set cr value
      cr <- 2*e - ty %% 2
      
      for (fdr in 1:nF) {
        for (y in 1:nY) {
        
        
        # biomass index
        index_B <- (ty - 1)*nE*nF*nY + (e - 1)*nF*nY + (fdr - 1)*nY + y 
        # yield index
        index_Y <- index_B + nTy*nE*nF*nY

        # calculate and store cumulative actual values 
        DF$Value[index_B] <- mean(Rel_biomass[1:y, cr, fdr])
        DF$Value[index_Y] <- mean(Rel_yield[1:y, cr, fdr])

        }
      }
    }
  }
  
  # make FDR and Estimate factor variables
  DF$FDR <- factor(DF$FDR)
  DF$Estimate <- factor(DF$Estimate, levels = c('Low', 'True', 'High'))
  DF$Type.Metric <- paste(DF$Type, DF$Metric, sep = ' ')
  DF$Type.Metric <- factor(DF$Type.Metric, 
                           levels = c('Static Biomass', 'Static Yield', 
                                      'Transient Biomass', 'Transient Yield'))
  
  # remove static biomass and yield for low and high estimates
  YEAR1 <- subset(DF, Year == years_to_plot[1] & 
                  (Estimate == 'True' | Type == 'Transient'))
  YEAR2 <- subset(DF, Year == years_to_plot[2] & 
                    (Estimate == 'True' | Type == 'Transient'))
  
  # panel labels
  panel.labels <- c(paste('(a) Biomass: year', years_to_plot[1]), 
                    paste('(b) Yield: year', years_to_plot[1]))
  names(panel.labels) <- c('Biomass', 'Yield')
  
  # plotting parameters
  og_colors <- rev(viridis(max(c(nF1, nF2)) + 1))
  if (s != 4) {
    new_colors <- og_colors[(nF2 - nF1 + 2):(nF2 + 1)]
  } else {
    new_colors <- og_colors[2:(nF2 + 1)]
  }
  
  # plot relative results
  thing1 <- ggplot(YEAR1, aes(x = Estimate, y = Value, color = FDR, 
                               shape = Type.Metric)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_point(position = position_jitter(w = 0.3, h = 0), size = 3, stroke = 1) +
    scale_shape_manual(values = c(16, 17, 1, 2), 
                       labels = c('Static Biomass', 'Static Yield', 
                                  'Transient Biomass', 'Transient Yield')) +
    scale_color_manual(values = new_colors, guide = 'none') +
    ylab('Cumulative value') +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    labs(shape = 'Type and Metric') +
    facet_grid(cols = vars(Metric), labeller = labeller(Metric = panel.labels))
  
  # panel labels
  panel.labels <- c(paste('(c) Biomass: year', years_to_plot[2]), 
                    paste('(d) Yield: year', years_to_plot[2]))  
  names(panel.labels) <- c('Biomass', 'Yield')
  
  # plot difference results
  thing2 <- ggplot(YEAR2, aes(x = Estimate, y = Value, color = FDR, 
                                shape = Type.Metric)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_point(position = position_jitter(w = 0.3, h = 0), size = 3, stroke = 1) +
    scale_shape_manual(values = c(16, 17, 1, 2), guide = 'none') +
    scale_color_manual(values = new_colors) +
    ylab('Cumulative Value') +
    xlab('Estimate of natural mortality (M)') +
    theme_bw() +
    labs(color = expression('D'[final])) +
    guides(color = guide_legend(order = 1)) +
    facet_grid(cols = vars(Metric), labeller = labeller(Metric = panel.labels))
  
  final_plot <- thing1 / thing2
  
  # save results to figures folder
  ggsave(final_plot, filename = paste(titles[s], 'M_', Names[s], '.png', sep = ''),
         path = 'C:/Users/vique/Box Sync/Quennessen_Thesis/MS thesis/publication manuscript/figures',
         width = png_width, height = png_height)
  
}
