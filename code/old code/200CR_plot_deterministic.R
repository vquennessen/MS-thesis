# relative biomass, yield, and effort

# plot relative biomass for static vs. transient DRs for each area

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(viridis)

###############################################################################
# CHECK THESE EVERY TIME
folder <- 'None'
cluster <- FALSE
png_width <- 4
png_height <- 6
alfa <- 0.25
all <- FALSE
###############################################################################

# species to compare
species_list <- c('CR_OR_2015')
Names <- c('Canary Rockfish')

# determine num_sims based on data folder
num_sims <- 3

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 200
if (all == TRUE) { 
  Final_DRs <- c(0.6, 0.7, 0.8, 0.9)
} else { 
  Final_DRs <- c(0.8, 0.9) }

types <- c('Static', 'Transient')
metrics <- c('Biomass', 'Yield', 'Effort', 'Density Ratio')

nY <- Time2 + 1
nS <- length(species_list)
nF <- length(Final_DRs)
nTy <- length(types)

DF <- data.frame(Type = rep(types, each = nF*nY), 
                 FDR = rep(Final_DRs, times = nTy, each = nY), 
                 Year = rep(0:Time2, times = nTy*nF), 
                 Value = rep(NA, nTy*nF*nY), 
                 Lower = rep(NA, nTy*nF*nY),
                 Upper = rep(NA, nTy*nF*nY))

# load biomass, yield, and effort files
folder <- ifelse(cluster == TRUE, 'Documents', 'Projects')

load(paste('~/', folder, '/MS-thesis/data/None/CR_OR_2015/200 year runs/', 
           '3_biomass.Rda', sep = ''))
load(paste('~/', folder, '/MS-thesis/data/None/CR_OR_2015/200 year runs/', 
           '3_yield.Rda', sep = ''))
load(paste('~/', folder, '/MS-thesis/data/None/CR_OR_2015/200 year runs/', 
           '3_effort.Rda', sep = ''))
load(paste('~/', folder, '/MS-thesis/data/None/CR_OR_2015/200 year runs/', 
           '3_DR.Rda', sep = ''))

##### relative biomass and median, upper, and lower limits  #####

# pull out sample sims as sums across all areas for particular years
B_sample   <- colSums(sims_biomass) 
Y_sample <- sims_yield
E_sample <- sims_effort

# initialize relative arrays
Rel_biomass <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
Rel_yield   <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
Rel_effort  <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))

# calculate relative arrays after reserve implementation
for (ty in 1:nTy) {
  for (fdr in 1:nF) {
    for (sim in 1:num_sims) {
      
      ind <- ifelse(all == TRUE, fdr, fdr + 2)
      
      Rel_biomass[, ty, fdr, sim] <- B_sample[, ty, ind, sim] / 
        B_sample[1, ty, ind, sim]
      Rel_yield[, ty, fdr, sim] <- Y_sample[, ty, ind, sim] / 
        Y_sample[1, ty, ind, sim]
      Rel_effort[, ty, fdr, sim] <- E_sample[Time1:(Time1 + Time2), 
                                             ty, ind, sim] / 
        E_sample[1, ty, ind, sim]
    }
  }
}

BIOMASS <- DF; YIELD <- DF; EFFORT <- DF; DR <- DF

##### fill in data frames with median and quantile values #####
for (ty in 1:nTy) {
  for (fdr in 1:nF) {
    for (y in 1:nY) {
      index <- (ty - 1)*nF*nY + (fdr - 1)*nY + y
      
      ind <- ifelse(all == TRUE, fdr, fdr + 2)
      
      # relative biomass
      BIOMASS$Value[index] <- median(Rel_biomass[y, ty, fdr, ])
      BIOMASS$Lower[index] <- quantile(Rel_biomass[y, ty, fdr, ], 0.25)
      BIOMASS$Upper[index] <- quantile(Rel_biomass[y, ty, fdr, ], 0.75)
      
      # relative yield
      YIELD$Value[index] <- median(Rel_yield[y, ty, fdr, ])
      YIELD$Lower[index] <- quantile(Rel_yield[y, ty, fdr, ], 0.25)
      YIELD$Upper[index] <- quantile(Rel_yield[y, ty, fdr, ], 0.75)
      
      # relative effort
      EFFORT$Value[index] <- median(Rel_effort[y, ty, fdr, ])
      EFFORT$Lower[index] <- quantile(Rel_effort[y, ty, fdr, ], 0.25)
      EFFORT$Upper[index] <- quantile(Rel_effort[y, ty, fdr, ], 0.75)
      
      DR$Value[index] <- median(sims_DR[y, ty, ind, ])
      DR$Lower[index] <- median(sims_DR[y, ty, ind, ])
      DR$Upper[index] <- median(sims_DR[y, ty, ind, ])
      
    }
  }
}

##### plotting parameters #####
jitter_height <- 0.006

if (all == TRUE) {
  og_colors <- viridis(5)
  new_colors <- og_colors[1:4]
} else {
  og_colors <- viridis(3)
  new_colors <- og_colors[1:2]
}


##### plot total biomass #####
biomass <- ggplot(data = BIOMASS, aes(x = Year, y = Value, 
                                      color = as.factor(FDR), 
                                      linetype = as.factor(Type))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR),
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(c(new_colors), alfa)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_color_manual(values = new_colors) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative biomass') +
  theme(axis.title.x = element_blank()) +
  labs(tag = 'a') +
  theme_bw() +
  theme(legend.position = 'none')

##### plot total yield #####
yield <- ggplot(data = YIELD, aes(x = Year, y = Value, color = as.factor(FDR), 
                                  linetype = as.factor(Type))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR),
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(c(new_colors), alfa)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_color_manual(values = new_colors) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative yield') +
  theme(axis.title.x = element_blank()) +
  labs(tag = 'b') +
  theme_bw() +
  theme(legend.position = 'none')


##### plot total effort #####
effort <- ggplot(data = EFFORT, aes(x = Year, y = Value, 
                                    color = as.factor(FDR), 
                                    linetype = as.factor(Type))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR),
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(c(new_colors), alfa)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height*12)) +
  scale_color_manual(values = new_colors) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative effort') +
  xlab('Years since reserve implemented') +
  labs(color = expression('D'[final]), linetype = 'Type', tag = 'c') +
  theme_bw() +
  theme(plot.margin = unit(c(0, 70, 0, 0), 'pt')) +
  theme(legend.position = c(1.24, 2)) + 
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))

##### plot density ratio #####
dr <- ggplot(data = DR, aes(x = Year, y = Value, color = as.factor(FDR), 
                            linetype = as.factor(Type))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR),
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(c(new_colors), alfa)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height*12)) +
  scale_color_manual(values = new_colors) +
  ylab('Density Ratio') +
  xlab('Years since reserve implemented') +
  labs(color = expression('D'[final]), linetype = 'Type', tag = 'c') +
  theme_bw() +
  theme(legend.position = 'none')

##### patch all the figures together #####
patch <- biomass / yield / effort / dr

thing <- ifelse(all == TRUE, 'all', 'conservative')

if (cluster == TRUE) {
  ggsave(patch,
         filename = paste('200CR_deterministic_', thing, '.png', sep = ''),
         path = paste('~/Documents/MS-thesis/figures/200 year runs/', folder, 
                      sep = ''),
         width = png_width, height = png_height)
  
} else {
  ggsave(patch, filename = paste('200CR_deterministic_', thing, '.png', 
                                 sep = ''),
         path = paste('C:/Users/Vic/Box/Quennessen_Thesis/MS thesis/', 
                      'publication manuscript/figures/', sep = ''),
         width = png_width, height = png_height)
}

