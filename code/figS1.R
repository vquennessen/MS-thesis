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

###############################################################################
# CHECK THESE EVERY TIME
folder <- 'None'
cluster <- FALSE

png_width <- 5.5
png_height <- 5.5

###############################################################################

# species to compare
species_list <- c('CR_OR_2015')
Names <- c('Canary Rockfish')

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 100
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:6)
types <- c('Static', 'Transient')
metrics <- c('Biomass', 'Yield', 'Effort')
MSY_metrics <- c('Biomass')
values <- c(0.4)

# dimensions
num_sims <- 1
nT <- Time2 + 1
nC <- length(Control_rules)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nTy <- length(types)
nM <- length(metrics)

base1 <- data.frame(Type = rep(types, each = nF1*nT), 
                    FDR = rep(Final_DRs1, times = 2, each = nT), 
                    Year = rep(0:Time2, times = 2*nF1), 
                    Value = rep(NA, 2*nF1*nT), 
                    Lower = rep(NA, 2*nF1*nT),
                    Upper = rep(NA, 2*nF1*nT))

for (s in 1:length(species_list)) {
  
  # load biomass, yield, and effort files
  if (cluster == TRUE) {
    load(paste('~/Documents/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_effort.Rda', sep = ''))
    
  } else {
    load(paste('~/Projects/MS-thesis/data/', folder, '/', 
               species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
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
  
  ##### fill in data frames with median and quantile values #####
  for (ty in 1:2) {
    for (fdr in 1:nF) {
      for (t in 1:nT) {
        index <- (ty - 1)*nF*nT + (fdr - 1)*nT + t
        
        # relative biomass
        BIOMASS$Value[index] <- median(Rel_biomass[t, ty, fdr])
        BIOMASS$Lower[index] <- quantile(Rel_biomass[t, ty, fdr], 0.25)
        BIOMASS$Upper[index] <- quantile(Rel_biomass[t, ty, fdr], 0.75)
        
        # relative yield
        YIELD$Value[index] <- median(Rel_yield[t, ty, fdr])
        YIELD$Lower[index] <- quantile(Rel_yield[t, ty, fdr], 0.25)
        YIELD$Upper[index] <- quantile(Rel_yield[t, ty, fdr], 0.75)
        
        # relative effort
        EFFORT$Value[index] <- median(Rel_effort[t, ty, fdr])
        EFFORT$Lower[index] <- quantile(Rel_effort[t, ty, fdr], 0.25)
        EFFORT$Upper[index] <- quantile(Rel_effort[t, ty, fdr], 0.75)
        
      }
    }
  }
  
  # put dataframes together, with new metric column
  BIOMASS$Metric <- 'Biomass'
  YIELD$Metric <- 'Yield'
  EFFORT$Metric <- 'Effort'
  DF <- rbind(BIOMASS, YIELD, EFFORT)
  DF$Metric <- factor(DF$Metric, levels = metrics)
  
  # calculate MSY
  source("code/calculate_MSY.R")
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
  
  ##### plotting parameters #####
  og_colors <- rev(viridis(max(c(nF1, nF2)) + 1))
  if (s != 4) {
    new_colors <- og_colors[(nF2 - nF1 + 2):(nF2 + 1)]
  } else {
    new_colors <- og_colors[2:(nF2 + 1)]
  }
  
  ##### new plot #####
  fig <- ggplot(data = DF, aes(x = Year, y = Value, color = as.factor(FDR), 
                               linetype = as.factor(Type))) +
    geom_hline(data = MSY_DF, aes(yintercept = Value, color = 'MSY'), 
               size = 0.75, linetype = 'twodash', color = 'black') +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    geom_line(size = 0.75) +
    scale_color_manual(values = new_colors) +
    facet_grid(Metric ~ Type, scales = 'free', switch = 'y') +
    ylab('Relative Value') +
    labs(color = expression('D'[final]), 
         linetype = 'Type of \n Control \n Rule') +
    theme_bw()
  
  
  # add panel tags (a) through (f)
  final_plot <- tag_facet(p = fig, 
                          hjust = -0.1, 
                          vjust = 13.25) +    
    theme(strip.text = element_text(), strip.background = element_rect())
  
  ggsave(final_plot, 
         filename = 'figS1_Canary Rockfish.png',
         path = 'C:/Users/Vic/Box Sync/Quennessen_Thesis/MS thesis/publication manuscript/figures',
         width = png_width, height = png_height)
  
}
