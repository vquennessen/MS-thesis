# relative biomass, yield, and effort

# plot relative biomass for static vs. transient DRs for each area

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(viridis)
library(egg)
library(cowplot)

###############################################################################
# CHECK THESE EVERY TIME
folder <- 'None'
cluster <- FALSE

png_width <- 5
png_height <- 5
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
# species_list <- c('CR_OR_2015_SSS', 'BR_OR_2015_SSS', 'LING_OW_2017_SSS', 
#                   'CAB_OR_2019_SSS')

Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')
titles <- c('old_figS1_', 'fig3_', 'figS2_', 'figS3_')
# titles <- c('SSS_old_figS1_', 'SSS_fig3_', 'SSS_figS2_', 'SSS_figS3_')


# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:6)
types <- c('Static', 'Transient')
metrics <- c('Biomass', 'Yield', 'Effort')
MSY_metrics <- c('Biomass', 'Biomass', 'SSB', 'Biomass')
values <- c(0.4, 0.4, 0.4, 0.4)

# dimensions
num_sims <- 1
nT <- Time2 + 1
nC <- length(Control_rules)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1
nTy <- length(types)
nM <- length(metrics)

# plotting things
MSY_x <- c(0.1675, 0.19, 0.18, 0.22)
MSY_y <- c(0.71, 0.72, 0.735, 0.715)
facets_v <- c(1.5, 3, 2, 11)

# initialise dataframes
base1 <- data.frame(Type = rep(types, each = nF1*nT), 
                    FDR = rep(Final_DRs1, times = nTy, each = nT), 
                    Year = rep(0:Time2, times = nTy*nF1), 
                    Value = rep(NA, nTy*nF1*nT)) 

base2 <- data.frame(Type = rep(types, each = nF2*nT), 
                    FDR = rep(Final_DRs2, times = nTy, each = nT), 
                    Year = rep(0:Time2, times = nTy*nF2), 
                    Value = rep(NA, nTy*nF2*nT))

for (s in 1:length(species_list)) {
  
  # load biomass, yield, and effort files
  if (cluster == TRUE) {
    load(paste('~/Documents/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_SSB.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_effort.Rda', sep = ''))
    
  } else {
    load(paste('~/Projects/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
    load(paste('~/Projects/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_SSB.Rda', sep = ''))
    load(paste('~/Projects/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
    load(paste('~/Projects/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_effort.Rda', sep = ''))
  }
  
  # set nF value for species 
  nF <- ifelse(s == 4, nF2, nF1)  

  ##### relative biomass and median, upper, and lower limits  #####
  
  # pull out sample sims as sums across all areas for particular years
  B_sample   <- colSums(sims_biomass) 
  SSB_sample   <- colSums(sims_SSB)
  Y_sample <- sims_yield
  E_sample <- sims_effort
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, nT*nC*nF), c(nT, nC, nF))
  Rel_yield   <- array(rep(0, nT*nC*nF), c(nT, nC, nF))
  Rel_effort  <- array(rep(0, nT*nC*nF), c(nT, nC, nF))
  
  # calculate relative arrays after reserve implementation
  for (cr in 1:nC) {
    for (fdr in 1:nF) {
      Rel_biomass[, cr, fdr] <- B_sample[1:nT, cr, fdr, 1] / 
        B_sample[1, cr, fdr, 1]
      Rel_yield[, cr, fdr] <- Y_sample[1:nT, cr, fdr, 1] / 
        Y_sample[1, cr, fdr, 1]
      Rel_effort[1:nT, cr, fdr] <- E_sample[Time1:(Time1 + Time2), cr, fdr, 1] / 
        E_sample[1, cr, fdr, 1]
    }
  }
  
  # data frame based on species number
  if (s == 4) { DF <- base2 } else { DF <- base1 }
  
  BIOMASS <- DF; YIELD <- DF; EFFORT <- DF
  
  # calculate MSY
  source("~/Projects/MS-thesis/code/calculate_MSY.R")
  MSY <- calculate_MSY(species_list[s], metric = MSY_metrics[s], 
                       value = values[s])
  
  MSY_biomass <- MSY[[2]]
  MSY_yield <- MSY[[4]]
  
  # calculate relative MSY values
  relative_MSY_biomass <- MSY_biomass / B_sample[1, 1, 1, 1]
  relative_MSY_yield <- MSY_yield / Y_sample[1, 1, 1, 1]
  
  MSY_values <- c(relative_MSY_biomass, relative_MSY_yield, NA)
  
  # MSY dataframe
  MSY_DF <- data.frame(Metric = rep(metrics, each = nTy), 
                       Types = rep(types, times = nM), 
                       Value = rep(MSY_values, each = nTy))
  
  MSY_DF$Metric <- factor(MSY_DF$Metric, levels = metrics)
  
  ##### fill in data frames with median and quantile values #####
  for (ty in 1:2) {
    for (fdr in 1:nF) {
      for (t in 1:nT) {
        index <- (ty - 1)*nF*nT + (fdr - 1)*nT + t
        
        # relative biomass
        BIOMASS$Value[index] <- median(Rel_biomass[t, ty, fdr])
        
        # relative yield
        YIELD$Value[index] <- median(Rel_yield[t, ty, fdr])
        
        # relative effort
        EFFORT$Value[index] <- median(Rel_effort[t, ty, fdr])
        
      }
    }
  }

  # put dataframes together, with new metric and MSY column
  BIOMASS$Metric <- 'Biomass'
  YIELD$Metric <- 'Yield'
  EFFORT$Metric <- 'Effort'
  DF <- rbind(BIOMASS, YIELD, EFFORT)
  DF$Metric <- factor(DF$Metric, levels = metrics)
  
  ##### plotting parameters #####
  og_colors <- rev(viridis(max(c(nF1, nF2)) + 1))
  if (s != 4) {
    new_colors <- og_colors[(nF2 - nF1 + 2):(nF2 + 1)]
  } else {
    new_colors <- og_colors[2:(nF2 + 1)]
  }
  
  ##### new plot #####
  fig1 <- ggplot(data = DF, aes(x = Year, y = Value, color = as.factor(FDR))) +
    geom_hline(data = MSY_DF, aes(yintercept = Value), 
               size = 0.75, linetype = 'twodash') +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black', 
               size = 0.5) +
    geom_line(size = 0.75) +
    scale_color_manual(values = c(new_colors)) +
    facet_grid(Metric ~ Type, scales = 'free', switch = 'y') +
    ylab('Relative Value') +
    labs(color = expression('D'[final])) +
    theme_bw()
  
  # add panel tags (a) through (f)
  fig2 <- tag_facet(p = fig1, hjust = -0.3, vjust = facets_v[s]) +
    theme(strip.text = element_text(), strip.background = element_rect())  
  
  final_plot <- ggdraw(fig2) + 
    draw_label(label = "MSY", x = MSY_x[s], y = MSY_y[s], size = 8.5)
  
  ggsave(final_plot, filename = paste(titles[s], Names[s], '.png', 
                                      sep = ''),
         path = 'C:/Users/Vic/Box Sync/Quennessen_Thesis/MS thesis/publication manuscript/figures',
         width = png_width, height = png_height)
  
}
