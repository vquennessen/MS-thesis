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
###############################################################################

# species to compare
species_list <- c('BR_OR_2015')
Names <- c('Black Rockfish')

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
num_sims <- 1
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
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/1_biomass.Rda', sep = ''))
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/1_yield.Rda', sep = ''))
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/1_effort.Rda', sep = ''))
load(paste('~/Projects/MS-thesis/data/None/BR_OR_2015/1_DR.Rda', sep = ''))

##### relative biomass and median, upper, and lower limits  #####

# pull out sample sims as sums across all areas for particular years
B_sample   <- colSums(sims_biomass) 
Y_sample <- sims_yield
E_sample <- sims_effort
DR_sample <- sims_DR

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

# set base data frame 
BIOMASS <- DF
YIELD   <- DF
EFFORT  <- DF
DR      <- DF

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
      
      # relative density ratio
      DR$Value[index] <- median(DR_sample[t, ty, fdr, 1])
      DR$Lower[index] <- quantile(DR_sample[t, ty, fdr, 1], 0.25)
      DR$Upper[index] <- quantile(DR_sample[t, ty, fdr, 1], 0.75)
      
    }
  }
}

##### plotting parameters #####
jitter_height <- 0
plot_color <- 'black'

##### single dfinal value plot #####
BIOMASS <- subset(BIOMASS, FDR == 0.9)
YIELD <- subset(YIELD, FDR == 0.9)
EFFORT <- subset(EFFORT, FDR == 0.9)
DR <- subset(DR, FDR == 0.9)

##### plot total biomass #####
biomass <- ggplot(data = BIOMASS, aes(x = Year, y = Value, 
                                      linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height), size = 1, 
            color = plot_color) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative biomass') +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = 'none') +
  annotate('text', x = 0, y = 0.957, label = 'a')

##### plot total yield #####
yield <- ggplot(data = YIELD, aes(x = Year, y = Value,
                                  linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height), size = 1, 
            color = plot_color) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative yield') +
  theme_bw() +
  theme(legend.position = 'none') +
  annotate('text', x = 0, y = 0.52, label = 'b') +
  theme(axis.title.x = element_blank())

##### plot total effort #####
effort <- ggplot(data = EFFORT, aes(x = Year, y = Value, 
                                    linetype = as.factor(Type))) +
  geom_line(position = position_jitter(w = 0, h = jitter_height), size = 1, 
            color = plot_color) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ylab('Relative effort') +
  theme_bw() +
  theme(legend.position = 'none') + 
  annotate('text', x = 0, y = 0.475, label = 'c') +
  theme(axis.title.x = element_blank())

# transient target density ratio values
years <- 0:20
M <- parameters(species_list[1])[[2]]
targets <- 1 - (1 - FDRs1[1])*(1 - exp(-1 * M * years))
transients <- data.frame(Year = years, Target = targets)

# add targets to DR dataframe
DR$Type2 <- c(rep('Static Actual', Time2 + 1), 
              rep('Transient Actual', Time2 + 1))
DR$Type3 <- 'Actual'
DR2 <- data.frame(Type = c(rep('Static', Time2 + 1), 
                           rep('Transient', Time2 + 1)), 
                  FDR = 0.9, 
                  Year = rep(0:20, times = 2), 
                  Value = c(rep(0.9, Time2 + 1), 
                            transients$Target), 
                  Lower = NA, Upper = NA, 
                  Type2 = c(rep('Static Target', Time2 + 1), 
                            rep('Transient Target', Time2 + 1)),
                  Type3 = 'Target')

DR3 <- rbind(DR, DR2)

##### plot total density ratio #####
dr <- ggplot(data = DR3, aes(x = Year, y = Value, linetype = Type2, 
                             color = Type2)) +
  geom_line(size = 0.8) +
  scale_color_manual(values = c('black', 'gray70', 'black', 'gray70')) +
  scale_linetype_manual(values = c(1, 1, 2, 2)) +
  ylab('Density ratio') +
  xlab('Years since reserve implemented') +
  scale_y_continuous(limits = c(0.84, 1.3)) +
  theme_bw() +
  theme(legend.position = c(0.5, 0.7)) +
  guides(color = guide_legend(nrow = 2)) +
  theme(legend.key.width = unit(1.25, 'cm')) +
  annotate('text', x = 0, y = 0.86, label = 'd') +
  labs(color = 'Control Rules', linetype = 'Control Rules')

##### patch all the figures together #####
patch <- biomass / yield / effort / dr

ggsave(patch, filename = paste('fig2_', Time2, '_years.png', sep = ''),
       path = 'C:/Users/Vic/Box Sync/Quennessen_Thesis/MS Thesis/publication manuscript/figures',
       width = png_width, height = png_height)
