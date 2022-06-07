# relative biomass, yield, and effort

# plot relative biomass for static vs. transient DRs for each area

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(egg)
library(viridis)
library(cowplot)

###############################################################################
# CHECK THESE EVERY TIME

folders <- c('Sampling', 'Recruitment', 'Both')
scenarios <- c('Sampling', 'Recruitment', 'Both')
versions <- c(1, 2)
max_FDRs <- c(0.8, 0.7, 0.8, 0.5)
FDRs1 <- c(0.6, 0.9)
FDRs2 <- c(0.4, 0.9)
cluster <- TRUE
png_width <- 7
png_height <- 6
alfa <- 0.25
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')
titles <- c('figS12_', 'figS11_', 'figS13_', 'figS14_',
            'figS8_', 'fig6_', 'figS9_', 'figS10_')

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

nSc <- length(scenarios)
nM  <- length(metrics)
nC  <- length(Control_rules)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFt <- 3*nF1 + nF2
nY  <- Time2 + 1
nTy <- length(types)

base1 <- data.frame(Species = rep(NA, nSc*nM*2*nF1*nY),
                    Scenario = rep(scenarios, each = nM*2*nF1*nY),
                    Metric = rep(metrics, times = nSc, each = 2*nF1*nY),
                    Type = rep(types, times = nSc*nM, each = nF1*nY), 
                    FDR = rep(Final_DRs1, times = nSc*nM*2, each = nY), 
                    Year = rep(0:Time2, times = nSc*nM*2*nF1), 
                    Value = rep(NA, nSc*nM*2*nF1*nY), 
                    Lower = rep(NA, nSc*nM*2*nF1*nY),
                    Upper = rep(NA, nSc*nM*2*nF1*nY), 
                    Source = rep('Stochastic', nSc*nM*2*nF1*nY))

base2 <- data.frame(Species = rep(NA, nSc*nM*2*nF2*nY),
                    Scenario = rep(scenarios, each = nM*2*nF2*nY),
                    Metric = rep(metrics, times = nSc, each = 2*nF2*nY),
                    Type = rep(types, times = nSc*nM, each = nF2*nY), 
                    FDR = rep(Final_DRs2, times = nSc*nM*2, each = nY), 
                    Year = rep(0:Time2, times = nSc*nM*2*nF2), 
                    Value = rep(NA, nSc*nM*2*nF2*nY), 
                    Lower = rep(NA, nSc*nM*2*nF2*nY),
                    Upper = rep(NA, nSc*nM*2*nF2*nY), 
                    Source = rep('Stochastic', nSc*nM*2*nF2*nY))

for (v in 1:2) {
  
  if (v == 1) { 
    FDRs <- list(c(0.6, 0.9), c(0.6, 0.9), c(0.6, 0.9), c(0.4, 0.9))
    
  } else if (v == 2) {
    FDRs <- c(0.8, 0.7, 0.8, 0.5)
  }
  
  for (s in 1:length(species_list)) {
    
    # data frame based on species number
    if (s == 4) { DF <- base2 } else { DF <- base1 }
    
    # add in species name
    DF$Species <- Names[s]
    
    for (f in 1:length(folders)) {
      
      num_sims <- 5000
      
      # load biomass, yield, and effort files
      load(paste('~/Documents/MS-thesis/data/', folders[f], '/', 
                 species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
      load(paste('~/Documents/MS-thesis/data/', folders[f], '/', 
                 species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
      load(paste('~/Documents/MS-thesis/data/', folders[f], '/', 
                 species_list[s], '/', num_sims, '_effort.Rda', sep = ''))
      
      # set nF value for species 
      nF <- ifelse(s == 4, nF2, nF1)  
      
      ##### relative biomass and median, upper, and lower limits  #####
      
      # pull out sample sims as sums across all areas for particular years
      B_sample   <- colSums(sims_biomass[, , 1:2, , ]) 
      Y_sample <- sims_yield[ , 1:2, , ]
      E_sample <- sims_effort[ , 1:2, , ]
      
      # initialize relative arrays
      Rel_biomass <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
      Rel_yield   <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
      Rel_effort  <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
      
      # calculate relative arrays after reserve implementation
      for (ty in 1:nTy) {
        for (fdr in 1:nF) {
          for (sim in 1:num_sims) {
            Rel_biomass[, ty, fdr, sim] <- B_sample[1:nY, ty, fdr, sim] / 
              B_sample[1, ty, fdr, sim]
            Rel_yield[, ty, fdr, sim] <- Y_sample[1:nY, ty, fdr, sim] / 
              Y_sample[1, ty, fdr, sim]
            Rel_effort[, ty, fdr, sim] <- E_sample[Time1:(Time1 + Time2), 
                                                   ty, fdr, sim] / 
              E_sample[1, ty, fdr, sim]
          }
        }
      }
      
      # vectors of relative values
      vectors <- list(Rel_biomass, Rel_yield, Rel_effort)
      
      ##### fill in data frames with median and quantile values #####
      for (m in 1:nM) {
        for (ty in 1:nTy) {
          for (fdr in 1:nF) {
            for (y in 1:nY) {
              
              index <- (f - 1)*nM*2*nF*nY + (m - 1)*2*nF*nY + (ty - 1)*nF*nY + 
                (fdr - 1)*nY + y
              # print(index)
              
              # set relative vector
              vector <- vectors[[m]]
              
              # relative biomass
              DF$Value[index] <- median(vector[y, ty, fdr, ])
              DF$Lower[index] <- quantile(vector[y, ty, fdr, ], 0.25)
              DF$Upper[index] <- quantile(vector[y, ty, fdr, ], 0.75)
              
            }
          }
        }
        
      }
      
    }
    
    # calculate MSY
    source("calculate_MSY.R")
    MSY <- calculate_MSY(species_list[s], metric = MSY_metrics[s], 
                         value = values[s])
    
    MSY_biomass <- MSY[[2]]
    MSY_yield <- MSY[[4]]
    
    # calculate relative MSY values
    relative_MSY_biomass <- MSY_biomass / B_sample[1, 1, 1, 1]
    relative_MSY_yield <- MSY_yield / Y_sample[1, 1, 1, 1]
    
    MSY_values <- c(relative_MSY_biomass, relative_MSY_yield, NA)
    
    # MSY dataframe
    MSY_DF <- data.frame(Metric = rep(metrics, each = nSc), 
                         Scenario = rep(scenarios, times = nM), 
                         Value = rep(MSY_values, each = nSc))
    
    # make columns match up with DF
    MSY_DF$Metric <- factor(MSY_DF$Metric, levels = metrics, 
                            labels = c('Relative Biomass', 
                                       'Relative Yield', 
                                       'Relative Effort'))
    MSY_DF$Scenario <- factor(MSY_DF$Scenario, levels = scenarios)
    
    
    ##### plotting parameters #####
    colors <- viridis(3)
    new_colors <- colors[1:2]
    
    # process DF
    DF$Scenario <- factor(DF$Scenario, levels = scenarios)
    DF$Metric <- factor(DF$Metric, levels = metrics, 
                        labels = c('Relative Biomass', 'Relative Yield', 
                                   'Relative Effort'))
    DF$FDR <- factor(DF$FDR, levels = Final_DRs2)
    DF$Type <- factor(DF$Type, levels = types)
    new_DF <- subset(DF, FDR %in% FDRs[[s]])
    
    if (v == 1) {
      
      fig1 <- ggplot(data = new_DF, aes(x = Year, y = Value, 
                                              linetype = Type, color = FDR)) +
        geom_hline(data = MSY_DF, aes(yintercept = Value), 
                   size = 0.75, linetype = 'twodash') +
        geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
        geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = FDR), 
                    color = NA, show.legend = FALSE) +
        scale_fill_manual(values = alpha(new_colors, alfa)) +
        geom_line() +
        scale_color_manual(values = new_colors) +
        theme_bw() +
        theme(axis.title.y = element_blank()) +
        facet_grid(Metric ~ Scenario, scales = 'free', switch = 'y') +
        xlab('Years since reserve implemented') +
        guides(color = guide_legend(order = 1), 
               linetype = guide_legend(order = 2)) +
        labs(color = expression('D'[final]), 
             linetype = 'Type of \n Control Rule')
      
      final_plot <- ggdraw(fig1) + 
        draw_label(label = "Reference", x = .9, y = .831, size = 11.5) + 
        draw_line(x = c(.86, .9), y = c(.785,.785), linetype = "twodash", 
                  size = .8) +
        draw_label(label = "MSY", x = .93, y = .785, size = 8.5) +
        draw_line(x = c(.86, .9), y = c(.745,.745), linetype = "dashed", 
                  size = .5) +
        draw_label(label = "No change", x = .954, y = .745, size = 8.5) 
      
    } else {
      
      ##### load deterministic values and merge into new dataframe #####
      
      # load deterministic data for correct FDR values
      source('~/Documents/MS-thesis/code/extract_deterministic_values.R')
      deterministic <- read.csv('deterministic_values.csv')
      new_deterministic <- subset(deterministic, 
                                  FDR == FDRs[s] & Species == Names[s])
      
      # combine stochastic and deterministic dataframes into one
      combined_DF <- rbind(new_DF, new_deterministic)
      
      # write combined_DF to .csv file
      write.csv(x = combined_DF, file = 'combinedDF1.csv')
      
      # process DF
      combined_DF$Scenario <- factor(combined_DF$Scenario, levels = scenarios)
      combined_DF$Metric <- factor(combined_DF$Metric)
      # , levels = metrics, 
      # labels = c('Relative Biomass', 
      #            'Relative Yield', 
      #            'Relative Effort'))
      
      # write combined_DF to .csv file
      write.csv(x = combined_DF, file = 'combinedDF2.csv')
      
      ##### plot panels ##########################################################
      fig1 <- ggplot(data = combined_DF, 
                           aes(x = Year, y = Value, color = Type, 
                               linetype = Source)) +
        geom_hline(data = MSY_DF, aes(yintercept = Value), 
                   size = 0.75, linetype = 'twodash') +
        geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
        geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Type), 
                    color = NA, show.legend = FALSE) +
        scale_fill_manual(values = alpha(new_colors, alfa)) +
        geom_line() +
        scale_color_manual(values = new_colors) +
        theme_bw() +
        theme(axis.title.y = element_blank()) +
        facet_grid(Metric ~ Scenario, scales = 'free', switch = 'y') +
        xlab('Years since reserve implemented') +
        guides(color = guide_legend(order = 1), 
               linetype = guide_legend(order = 2)) +
        labs(color = 'Type of \n Control Rule', linetype = 'Source')
      
      final_plot <- ggdraw(fig1) + 
        draw_label(label = "Reference", x = .884, y = .82, size = 11.5) + 
        draw_line(x = c(.84, .875), y = c(.775,.775), linetype = "twodash", 
                  size = .8) +
        draw_label(label = "MSY", x = .905, y = .775, size = 8.5) +
        draw_line(x = c(.84, .875), y = c(.735,.735), linetype = "dashed", 
                  size = .5) +
        draw_label(label = "No change", x = .928, y = .735, size = 8.5) 
      
    }
    
    final_plot <- tag_facet(final_plot) +
      theme(strip.text = element_text(), strip.background = element_rect())
    
    # figure titles index
    fti <- (v - 1)*length(species_list) + s
    
    # save plot
    ggsave(final_plot,
           filename = paste(titles[fti], Names[s], '.png', sep = ''),
           path = '~/Documents/MS-thesis/figures/',
           width = png_width, height = png_height)
  }
  
}
