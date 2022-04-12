# relative biomass, yield, and effort

# plot relative biomass for static vs. transient DRs for each area

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
FDRs1 <- c(0.9)
cluster <- FALSE
png_width <- 4
png_height <- 8
y1 <- 0.25
y2 <- 1.5
y1.1 <- 0.95
y2.1 <- 2
alfa <- 0.25
###############################################################################

# species to compare
species_list <- c('BR_OR_2015')
Names <- c('Black Rockfish')

# determine num_sims based on data folder
num_sims <- 3

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:6)
types <- c('Static', 'Transient')
metrics <- c('Biomass', 'Yield', 'Effort', 'Density.Ratio')

# dimensions
sample_size = num_sims

nT <- Time2 + 1
nC <- length(Control_rules)
nS <- length(species_list)
nF <- length(Final_DRs1)

DF <- data.frame(Type = rep(types, each = nF*nT), 
                   FDR = rep(Final_DRs1, times = 2, each = nT), 
                   Year = rep(0:Time2, times = 2*nF), 
                   Value = rep(NA, 2*nF*nT), 
                   Lower = rep(NA, 2*nF*nT),
                   Upper = rep(NA, 2*nF*nT))

# load biomass, yield, and effort files
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/3_biomass.Rda', sep = ''))
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/3_yield.Rda', sep = ''))
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/3_effort.Rda', sep = ''))
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/3_DR.Rda', sep = ''))
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/3_abundance.Rda', sep = ''))

##### relative biomass and median, upper, and lower limits  #####

# pull out sample sims as sums across all areas for particular years
B_sample  <- colSums(sims_biomass) 
Y_sample  <- sims_yield
E_sample  <- sims_effort
DR_sample <- sims_DR
A_sample  <- colSums(sims_abundance)

# initialize relative arrays
Rel_biomass <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
Rel_yield   <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
Rel_effort  <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
Rel_abundance  <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))

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
        E_sample[Time1, cr, fdr, sim]
      Rel_abundance[, cr, fdr, sim] <- A_sample[, cr, fdr, sim] / 
        A_sample[1, cr, fdr, sim]
    }
  }
}

# set base data frame 
BIOMASS   <- DF
YIELD     <- DF
EFFORT    <- DF
DR        <- DF
ABUNDANCE <- DF

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
      
      # relative density ratio
      DR$Value[index] <- median(DR_sample[t, ty, fdr, ])
      DR$Lower[index] <- quantile(DR_sample[t, ty, fdr, ], 0.25)
      DR$Upper[index] <- quantile(DR_sample[t, ty, fdr, ], 0.75)
      
      # relative abundance
      ABUNDANCE$Value[index] <- median(Rel_abundance[t, ty, fdr, ])
      ABUNDANCE$Lower[index] <- quantile(Rel_abundance[t, ty, fdr, ], 0.25)
      ABUNDANCE$Upper[index] <- quantile(Rel_abundance[t, ty, fdr, ], 0.75)
      
    }
  }
}

##### plotting parameters #####
jitter_height <- 0

##### single dfinal value plot #####
BIOMASS <- subset(BIOMASS, FDR == 0.9)
YIELD <- subset(YIELD, FDR == 0.9)
EFFORT <- subset(EFFORT, FDR == 0.9)
DR <- subset(DR, FDR == 0.9)
ABUNDANCE <- subset(ABUNDANCE, FDR == 0.9)

##### plot total abundance #####
abundance <- ggplot(data = ABUNDANCE, aes(x = Year, y = Value, 
                                          linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height), size = 1) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative \n abundance') +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = 'none') +
  annotate('text', x = 0, y = 1.21, label = 'a')

##### plot total biomass #####
biomass <- ggplot(data = BIOMASS, aes(x = Year, y = Value, 
                                      linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height), size = 1) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative \n biomass') +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = 'none') +
  annotate('text', x = 0, y = 1.21, label = 'b')

##### plot total yield #####
yield <- ggplot(data = YIELD, aes(x = Year, y = Value,
                                  linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height), size = 1) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative yield') +
  theme_bw() +
  theme(legend.position = 'none') +
  annotate('text', x = 0, y = 1.15, label = 'c') +
  theme(axis.title.x = element_blank())

##### plot total effort #####
effort <- ggplot(data = EFFORT, aes(x = Year, y = Value, 
                                    linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height), size = 1) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative effort') +
  theme_bw() +
  theme(legend.position = 'none') + 
  annotate('text', x = 0, y = 1.9, label = 'd') +
  theme(axis.title.x = element_blank())

# transient target density ratio values
years <- 0:20
M <- parameters(species_list[1])[[2]]
targets <- 1 - (1 - FDRs1[1])*(1 - exp(-1 * M * years))
transients <- data.frame(Year = years, Target = targets)

##### plot total density ratio #####
dr <- ggplot(data = DR, aes(x = Year, y = Value, linetype = as.factor(Type))) +
  geom_hline(yintercept = 0.9, size = 0.75, color = 'gray60') +
  geom_line(data = transients, aes(x = Year, y = Target), linetype = 'dashed', 
            color = 'gray60', size = 0.75) + 
  geom_line(position = position_jitter(w = 0, h = jitter_height), size = 1) +
  ylab('Density ratio') +
  xlab('Years since reserve implemented') +
  labs(linetype = 'Type of Control Rule') +
  theme_bw() +
  # theme(plot.margin = unit(c(0, 70, 0, 0), 'pt')) +
  theme(legend.position = 'bottom') + # c(1.22, 2.28)) + 
  annotate('text', x = 0, y = 1.07, label = 'e')

##### patch all the figures together #####
patch <- abundance / biomass / yield / effort / dr 

ggsave(patch, filename = paste(Names[1], '_abundance.png', sep = ''),
       path = 'C:/Users/Vic/Box/Quennessen_Thesis/MS Thesis/publication manuscript/figures',
       width = png_width, height = png_height)
