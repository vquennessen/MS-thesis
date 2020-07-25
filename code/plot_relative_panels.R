# relative biomass, yield, and effort

# plot relative biomass for static vs. transient DRs for each area

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
folder <- 'None'
FDRs1 <- c(0.6, 0.9)
FDRs2 <- c(0.4, 0.9)
cluster <- FALSE
png_width <- 4
png_height <- 6
y1 <- 0.25
y2 <- 1.5
y1.1 <- 0.95
y2.1 <- 2
alfa <- 0.25
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
types <- c('Static', 'Transient')
metrics <- c('Biomass', 'Yield', 'Effort')

# dimensions
sample_size = num_sims

nT <- Time2 + 1
nC <- length(Control_rules)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1

base1 <- data.frame(Type = rep(types, each = nF1*nT), 
                    FDR = rep(Final_DRs1, times = 2, each = nT), 
                    Year = rep(0:Time2, times = 2*nF1), 
                    Value = rep(NA, 2*nF1*nT), 
                    Lower = rep(NA, 2*nF1*nT),
                    Upper = rep(NA, 2*nF1*nT))

base2 <- data.frame(Type = rep(types, each = nF2*nT), 
                    FDR = rep(Final_DRs2, times = 2, each = nT), 
                    Year = rep(0:Time2, times = 2*nF2), 
                    Value = rep(NA, 2*nF2*nT), 
                    Lower = rep(NA, 2*nF2*nT), 
                    Upper = rep(NA, 2*nF2*nT))

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
  
  # data frame based on species number
  if (s == 4) { DF <- base2 } else { DF <- base1 }
  
  BIOMASS <- DF; YIELD <- DF; EFFORT <- DF
  
  ##### fill in data frames with median and quantile values #####
  for (ty in 1:2) {
    for (fdr in 1:nF) {
      for (t in 1:nT) {
        index <- (ty - 1)*nF*nT + (fdr - 1)*nT + t
        
        # relative biomass
        BIOMASS$Value[index] <- median(Rel_biomass[t, ty, fdr, ])
        BIOMASS$Lower[index] <- quantile(Rel_biomass[t, ty, fdr, ], 0.25)
        BIOMASS$Upper[index] <- quantile(Rel_biomass[t, ty, fdr, ], 0.75)
        
        # relative yield
        YIELD$Value[index] <- median(Rel_yield[t, ty, fdr, ])
        YIELD$Lower[index] <- quantile(Rel_yield[t, ty, fdr, ], 0.25)
        YIELD$Upper[index] <- quantile(Rel_yield[t, ty, fdr, ], 0.75)
        
        # relative effort
        EFFORT$Value[index] <- median(Rel_effort[t, ty, fdr, ])
        EFFORT$Lower[index] <- quantile(Rel_effort[t, ty, fdr, ], 0.25)
        EFFORT$Upper[index] <- quantile(Rel_effort[t, ty, fdr, ], 0.75)
        
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
  
  # set FDRs for minimum and maximum biomass and yields and extract indices to 
  # get colors
  if (s == 4) { 
    ind <- which(Final_DRs2 %in% FDRs2)
  } else { 
    ind <- which(Final_DRs1 %in% FDRs1)
  }
  
  # pull colors out for scenarios that are not 'None'
  if (folder != 'None') {new_colors <- colors[ind]} else {new_colors <- colors}

  ##### plot total biomass #####
  biomass <- ggplot(data = subset(BIOMASS, FDR %in% FDRs), 
                    aes(x = Year, y = Value, color = as.factor(FDR), 
                        linetype = as.factor(Type))) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR),
                    colour = NA), show.legend = FALSE) +
    scale_fill_manual(values = alpha(c(new_colors), alfa)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_color_manual(values = new_colors) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    ylab('Relative Biomass') +
    theme(axis.title.x = element_blank()) +
    theme(legend.position = 'none')
  
  ##### plot total yield #####
  yield <- ggplot(data = subset(YIELD, FDR %in% FDRs), 
                  aes(x = Year, y = Value, color = as.factor(FDR), 
                      linetype = as.factor(Type))) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR),
                    colour = NA), show.legend = FALSE) +
    scale_fill_manual(values = alpha(c(new_colors), alfa)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_color_manual(values = new_colors) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    ylab('Relative Yield') +
    theme(legend.position = 'none') +
    theme(axis.title.x = element_blank())
  
  ##### plot total effort #####
  effort <- ggplot(data = subset(EFFORT, FDR %in% FDRs), 
                   aes(x = Year, y = Value, color = as.factor(FDR), 
                       linetype = as.factor(Type))) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR),
                    colour = NA), show.legend = FALSE) +
    scale_fill_manual(values = alpha(c(new_colors), alfa)) +
    geom_line(position = position_jitter(w = 0, h = jitter_height)) +
    scale_color_manual(values = new_colors) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
    ylab('Relative Effort') +
    labs(color = 'FDR', linetype = 'Type') +
    theme(plot.margin = unit(c(0, 60, 0, 0), 'pt')) +
    theme(legend.position = c(1.17, 1.5)) + 
    guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))
  
  ##### patch all the figures together #####
  patch <- biomass / yield / effort
  thing <- patch + plot_annotation(
    title = paste(Names[s], ': Relative Biomass, Yield, and Effort', sep = '')) &
    theme(text = element_text(size = 9))
  
  if (cluster == TRUE) {
    ggsave(thing, 
           filename = paste(Names[s], '_relative.png', sep = ''),
           path = paste('~/Documents/MS-thesis/figures/', folder, sep = ''),
           width = png_width, height = png_height)
    
  } else {
    ggsave(thing, filename = paste(Names[s], '_relative.png', sep = ''),
           path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                        folder, '/relative', sep = ''),
           width = png_width, height = png_height)
  }
  
}
