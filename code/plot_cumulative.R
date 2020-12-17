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
Years <- 1:20
png_width <- 8
png_height <- 6
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')

# determine num_sims based on data folder
num_sims <- 3

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = 1:6

# dimensions
sample_size = num_sims
PD = 0.25
types <- c('Static', 'Transient')
metrics <- c('Biomass', 'Yield')

nTy <- length(types)
nY  <- length(Years)
nM  <- length(metrics)
nC  <- length(Control_rules)
nS  <- length(species_list)
s1  <- 3; s2 <- 1
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFt <- 3*nF1 + nF2

DF1 <- data.frame(Species = c(rep(Names[1], each = nM*nTy*nF1*nY), 
                              rep(Names[2], each = nM*nTy*nF1*nY), 
                              rep(Names[3], each = nM*nTy*nF1*nY), 
                              rep(Names[4], each = nM*nTy*nF2*nY)),
                  Metric = c(rep(rep(metrics, each = nTy*nF1*nY), 3), 
                             rep(metrics, each = nTy*nF2*nY)),
                  Type = c(rep(rep(types, times = nM, each = nF1*nY), 3),
                           rep(types, times = nM, each = nF2*nY)),
                  FDR = c(rep(rep(Final_DRs1, times = nM*nTy, each = nY), 3),
                          rep(Final_DRs2, times = nM*nTy, each = nY)),
                  Year = rep(Years, times = nS*nM*nTy*nFt),
                  Value = rep(NA, nS*nM*nTy*nFt*nY))

for (s in 1:length(Names)) {
  
  # load biomass and yield files
  load(paste('~/Projects/MS-thesis/data/None/', 
             species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
  load(paste('~/Projects/MS-thesis/data/None/', 
             species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
  
  # set nF value for species 
  nF <- ifelse(s == 4, nF2, nF1)  
  
  ##### relative biomass and median, upper, and lower limits  #####
  
  # pull out sample sims as sums across all areas for particular years, but 
  # only for correct M value
  B_sample   <- colSums(sims_biomass[, , c(1, 4), , ]) 
  Y_sample <- sims_yield[, c(1, 4), , ]
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
  Rel_yield <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
  
  # calculate relative biomass and year arrays after reserve implementation
  for (ty in 1:nTy) {
    for (fdr in 1:nF) {
      for (y in 1:nY) {
        for (sim in 1:num_sims) {
          Rel_biomass[y, ty, fdr, sim] <- B_sample[y + 1, ty, fdr, sim] / 
            B_sample[1, ty, fdr, sim]
          Rel_yield[y, ty, fdr, sim] <- Y_sample[y + 1, ty, fdr, sim] / 
            Y_sample[1, ty, fdr, sim]
        }
      }
    }
  }
  
  ##### fill in data frame with cumulative values #####
  for (ty in 1:nTy) {
    for (fdr in 1:nF) {
      for (y in 1:nY) {
        
        # biomass
        
        # yield
      
      }
    }
  }
  
  # DF based on species
  if (s == 4) { diff_DF <- diff2 } else { diff_DF <- diff1 }
  
  # make FDR and Estimate factor variables
  DF$FDR <- factor(DF$FDR)
  DF$Type.Metric <- paste(DF$Type, DF$Metric, sep = ' ')
  DF$Type.Metric <- factor(DF$Type.Metric, 
                           levels = c('Static Biomass', 'Static Yield', 
                                      'Transient Biomass', 'Transient Yield'))
  
  # remove static biomass and yield for low and high estimates
  new_DF <- subset(DF, Estimate == 'True' | Type == 'Transient')
  
  # panel labels
  panel.labels <- c('(a) Biomass', '(b) Yield')
  names(panel.labels) <- c('Biomass', 'Yield')
  
  # plotting parameters
  og_colors <- rev(viridis(max(c(nF1, nF2)) + 1))
  if (s != 4) {
    new_colors <- og_colors[(nF2 - nF1 + 2):(nF2 + 1)]
  } else {
    new_colors <- og_colors[2:(nF2 + 1)]
  }
  
  # plot relative results
  thing1 <- ggplot(new_DF, aes(x = Estimate, y = Value, color = FDR, 
                               shape = Type.Metric)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_point(position = position_jitter(w = 0.25, h = 0), size = 3, stroke = 1) +
    scale_shape_manual(values = c(16, 17, 1, 2), 
                       labels = c('Static', 'Static', 'Transient', 'Transient')) +
    scale_color_manual(values = new_colors, guide = FALSE) +
    ylab('Median relative value') +
    theme_bw() +
    theme(axis.title.x = element_blank()) +
    labs(shape = 'Type', color = expression('D'[final])) +
    facet_grid(cols = vars(Metric), labeller = labeller(Metric = panel.labels))
  
  # make FDR and Estimate factor variables
  diff_DF$FDR <- factor(diff_DF$FDR)
  diff_DF$Estimate <- factor(diff_DF$Estimate, 
                             levels = c('Low', 'True', 'High'))
  
  # panel labels
  panel.labels <- c('(c) Biomass', '(d) Yield')
  names(panel.labels) <- c('Biomass', 'Yield')
  
  # plot difference results
  thing2 <- ggplot(diff_DF, aes(x = Estimate, y = Value, color = FDR, 
                                shape = Metric)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(position = position_jitter(w = 0.3, h = 0), size = 3, stroke = 1) +
    scale_shape_manual(values = c(1, 2), guide = FALSE) +
    scale_color_manual(values = new_colors) +
    ylab('Median difference') +
    xlab('Estimate of natural mortality (M)') +
    theme_bw() +
    labs(shape = 'Difference \n Metric', color = expression('D'[final])) +
    guides(color = guide_legend(order = 1)) +
    facet_grid(cols = vars(Metric), labeller = labeller(Metric = panel.labels))
  
  final_plot <- thing1 / thing2
  
  # save results to figures folder
  if (cluster == TRUE) {
    ggsave(final_plot, filename = paste(Names[s], '_year', Years[y], 
                                        '_M.png', sep = ''),
           path = paste('~/Documents/MS-thesis/figures/', folder, sep = ''),
           width = png_width, height = png_height)
    
  } else {
    ggsave(final_plot, filename = paste('M_', Names[s], '_year', Years[y], '.png', sep = ''),
           path = 'C:/Users/Vic/Box/Quennessen_Thesis/publication manuscript/viridis figures/',
           width = png_width, height = png_height)
  }
  
}