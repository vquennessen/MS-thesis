# relative biomass, yield, and effort

# plot relative biomass for static vs. transient DRs for each area

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
# install.packages('egg')
library(egg)

###############################################################################
# CHECK THESE EVERY TIME
folders <- c('Sampling', 'Both')
max_FDRs <- c(0.8, 0.7, 0.8, 0.5)
cluster <- TRUE
png_width <- 7
png_height <- 6
alfa <- 0.25
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
metrics <- c('Biomass', 'Yield', 'Effort')

nT <- Time2 + 1
nM <- length(metrics)
nC <- length(Control_rules)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)

for (s in 1:length(species_list)) {
  
  for (f in 1:length(folders)) {
    
    scenarios <- paste(folders[f], c('', '.HighT'), sep = '')
    nS <- length(scenarios)
    
    base1 <- data.frame(Scenario = rep(scenarios, each = nM*2*nF1*nT),
                        Metric = rep(metrics, times = nS, each = 2*nF1*nT),
                        Type = rep(types, times = nS*nM, each = nF1*nT), 
                        FDR = rep(Final_DRs1, times = nS*nM*2, each = nT), 
                        Year = rep(0:Time2, times = nS*nM*2*nF1), 
                        Value = rep(NA, nS*nM*2*nF1*nT), 
                        Lower = rep(NA, nS*nM*2*nF1*nT),
                        Upper = rep(NA, nS*nM*2*nF1*nT))
    
    base2 <- data.frame(Scenario = rep(scenarios, each = nM*2*nF2*nT),
                        Metric = rep(metrics, times = nS, each = 2*nF2*nT),
                        Type = rep(types, times = nS*nM, each = nF2*nT), 
                        FDR = rep(Final_DRs2, times = nS*nM*2, each = nT), 
                        Year = rep(0:Time2, times = nS*nM*2*nF2), 
                        Value = rep(NA, nS*nM*2*nF2*nT), 
                        Lower = rep(NA, nS*nM*2*nF2*nT),
                        Upper = rep(NA, nS*nM*2*nF2*nT))
    
    # data frame based on species number
    if (s == 4) { DF <- base2 } else { DF <- base1 }
    
    for (scen in 1:nS) {
      
      # determine num_sims based on data folder
      num_sims <- ifelse(scenarios[scen] == 'Both', 6193, 5000)
      
      # load biomass, yield, and effort files
      load(paste('~/Documents/MS-thesis/data/', scenarios[scen], '/', 
                 species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
      load(paste('~/Documents/MS-thesis/data/', scenarios[scen], '/', 
                 species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
      load(paste('~/Documents/MS-thesis/data/', scenarios[scen], '/', 
                 species_list[s], '/', num_sims, '_effort.Rda', sep = ''))
      
      # set nF value for species 
      nF <- ifelse(s == 4, nF2, nF1)  
      
      ##### relative biomass and median, upper, and lower limits  #####
      
      # pull out sample sims as sums across all areas for particular years
      B_sample   <- colSums(sims_biomass) 
      Y_sample <- sims_yield
      E_sample <- sims_effort
      
      # initialize relative arrays
      Rel_biomass <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
      Rel_yield   <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
      Rel_effort  <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
      
      # calculate relative arrays after reserve implementation
      for (cr in 1:nC) {
        for (fdr in 1:nF) {
          for (sim in 1:num_sims) {
            Rel_biomass[, cr, fdr, sim] <- B_sample[, cr, fdr, sim] / 
              B_sample[1, cr, fdr, sim]
            Rel_yield[, cr, fdr, sim] <- Y_sample[, cr, fdr, sim] / 
              Y_sample[1, cr, fdr, sim]
            Rel_effort[, cr, fdr, sim] <- E_sample[Time1:(Time1 + Time2), 
                                                   cr, fdr, sim] / 
              E_sample[1, cr, fdr, sim]
          }
        }
      }
      
      # vectors of relative values
      vectors <- list(Rel_biomass, Rel_yield, Rel_effort)
      
      ##### fill in data frames with median and quantile values #####
      for (m in 1:nM) {
        for (ty in 1:2) {
          for (fdr in 1:nF) {
            for (t in 1:nT) {
              
              index <- (scen - 1)*nM*2*nF*nT + (m - 1)*2*nF*nT + (ty - 1)*nF*nT + 
                (fdr - 1)*nT + t
              # print(index)
              
              # set relative vector
              vector <- vectors[[m]]
              
              # relative biomass
              DF$Value[index] <- median(vector[t, ty, fdr, ])
              DF$Lower[index] <- quantile(vector[t, ty, fdr, ], 0.25)
              DF$Upper[index] <- quantile(vector[t, ty, fdr, ], 0.75)
              
            }
          }
        }
        
      }
      
    }
    
    ##### plotting parameters #####
    jitter_height <- 0
    if (s != 4) {
      colors <- c("#00BA38", "#00BFC4", "#619CFF", "#F564E3")
    } else {
      colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
    }
    
    # set FDRs for maximum difference in biomass and/or yield for Both
    if (f == 2) {
      FDRs <- max_FDRs[s]
    } else if (f == 1 & s != 4) {
      FDRs <- Final_DRs1
    } else { FDRs <- Final_DRs2 }
    
  
    if (s != 4) { 
      ind <- which(Final_DRs1 %in% FDRs)
    } else if (s == 4) {
      ind <- which(Final_DRs2 %in% FDRs)
    }
    
    # pull colors out for scenarios that are not 'None'
    if (folders[f] != 'None') { new_colors <- colors[ind]
    } else { new_colors <- colors }
    
    # process DF
    levels(DF$Scenario) <- c('24 Transects', '36 Transects')
    DF$Metric <- factor(DF$Metric, levels = metrics)
    DF$FDR <- as.factor(DF$FDR)
    DF$Type <- as.factor(DF$Type)
    new_DF <- subset(DF, FDR %in% FDRs)
    
    ##### plot panels #####
    final_plot <- ggplot(data = new_DF, aes(x = Year, y = Value, color = FDR, 
                                            linetype = Type)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR),
                      colour = NA), show.legend = FALSE) +
      scale_fill_manual(values = alpha(c(new_colors), alfa)) +
      geom_line(position = position_jitter(w = 0, h = jitter_height)) +
      scale_color_manual(values = new_colors) +
      geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
      facet_grid(Metric ~ Scenario, scales = 'free') +
      labs(color = expression('D'[final]), linetype = 'Type')
    
    final_plot <- tag_facet(final_plot)
    final_plot <- final_plot +
      theme(strip.text = element_text(), strip.background = element_rect())
      
      
    # save plot
    ggsave(final_plot,
           filename = paste(Names[s], '_composite_relative.png', sep = ''),
           path = paste('~/Documents/MS-thesis/figures/', folders[f], '/', 
                        sep = ''),
           width = png_width, height = png_height)
    
  }
  
}
