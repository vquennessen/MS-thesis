# plot relative biomass compilation of all scenarios

# load any necessary libraries
# library(plyr)
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
png_width <- 7
png_height <- 7
y1 <- 0.25
y2 <- 1.5
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
species_names <- c('Canary \n Rockfish \n', 'Black \n Rockfish \n', 'Lingcod', 
                   'Cabezon')
nS <- length(species_list)

TimeT <- 100
years <- 0:TimeT
nT <- length(years)

types <- c('Static', 'Transient')

# estimates <- c('Low', 'True', 'High')
# nE <- length(estimates)

FDR <- 0.6

# DF <- data.frame(Species = rep(species_names, each = 2*nE*nT), 
#                  Type = rep(types, times = nS, each = nE*nT), 
#                  Estimate = rep(estimates, times = nS*2, each = nT),
#                  Time = rep(years, times = nS*2*nE), 
#                  Target = rep(0, nS*2*nE*nT))

DF <- data.frame(Species = rep(species_names, each = 2*nT), 
                 Type = rep(types, times = nS, each = nT), 
                 Time = rep(years, times = nS*2), 
                 Target = rep(0, nS*2*nT))

for (s in 1:nS) {
  
  # Natural mortality values
  par <- parameters(species_list[s])
  M <- par[[2]]
  # Nat_mortalities <- c(M - 0.05, M, M + 0.05)
  
  # static control rules
  static_start <- (s - 1)*2*nT + 1
  static_start
  static_end <- static_start + nT - 1
  static_end
  DF$Target[static_start:static_end] <- rep(FDR, nT)
  
  # for (e in 1:nE) {
  # 
  #   # natural mortality value
  #   NM <- Nat_mortalities[e]
  # 
  #   # transient control rules
  #   transient_start <- static_end + (e - 1)*nT  + 1
  #   transient_stop <- transient_start + nT - 1
  #   targets <- 1 - (1 - FDR)*(1 - exp(-1 * NM * years))
  #   DF$Target[transient_start:transient_stop] <- targets
  # 
  # }
    
    # transient control rules
    transient_start <- static_end + 1
    transient_start
    transient_stop <- transient_start + nT - 1
    transient_stop
    targets <- 1 - (1 - FDR)*(1 - exp(-1 * M * years))
    DF$Target[transient_start:transient_stop] <- targets
  
}

# make DF$Species a factor variable
DF$Species <- factor(DF$Species, levels = species_names)
# DF$Estimate <- factor(DF$Estimate, levels = estimates)

# dfinal text
dfinal <- expression('D'[final])

# plot it
# example <- ggplot(data = DF, aes(x = Time, y = Target, linetype = Estimate, 
#                                  size = Type)) +
#   geom_line() +
#   scale_size_manual(values = c(2, 1)) +
#   ylab('Target Density Ratio \n Dfinal value: 0.6') +
#   xlab('Years Since Reserve Implementation') +
#   scale_x_continuous(breaks = c(0, 10, 20)) +
#   ylim(0.5, 1) +
#   facet_grid(~ Species) +
#   labs(size = 'Control Rule \n Type', 
#        linetype = 'Estimate of \n Natural \n Mortality')

# only transient lines
DF_transient <- subset(DF, Type == 'Transient')

example <- ggplot(data = DF_transient, aes(x = Time, y = Target,
                                 linetype = Species)) +
  geom_line(size = 2) +
  xlab('Years Since Reserve Implementation') +
  scale_linetype_manual(values = c('solid', 'dotdash', 'dotted', 'dashed')) +
  # scale_color_manual(values = c('grey50', 'black')) +
  # scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_classic() +
  scale_size(range = 1, guide = 'none') +
  labs(linetype = 'Species', color = 'Type of \n Control Rule') +
  theme(legend.key.width = unit(1, 'cm')) +
  guides(linetype = guide_legend(override.aes = list(size = 1))) +
  theme(axis.title.y = element_blank())

ggsave(example, filename = 'fig1a_portion_transients_example.png', 
       path = 'C:/Users/Vic/Box Sync/Quennessen_Thesis/MS Thesis/publication manuscript/figures/',
       width = 20, height = 10)
