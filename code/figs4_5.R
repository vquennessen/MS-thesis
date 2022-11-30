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
png_width <- 9
png_height <- 6
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
# species_list <- c('CR_OR_2015_SSS', 'BR_OR_2015_SSS', 'LING_OW_2017_SSS', 
#                   'CAB_OR_2019_SSS')

Names <- c('Canary \n Rockfish', 'Black \n Rockfish', 'Lingcod', 'Cabezon')
years_to_plot <- c(10, 20)

# determine num_sims based on data folder
num_sims <- 1

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

DF <- data.frame(Species = c(rep(Names[1], each = nM*nTy*nF1*nY), 
                             rep(Names[2], each = nM*nTy*nF1*nY), 
                             rep(Names[3], each = nM*nTy*nF1*nY), 
                             rep(Names[4], each = nM*nTy*nF2*nY)),
                 Metric = c(rep(rep(metrics, each = nTy*nF1*nY), 3), 
                            rep(metrics, each = nTy*nF2*nY)),
                 Type = c(rep(rep(types, times = nM, each = nF1*nY), 3),
                          rep(types, times = nM, each = nF2*nY)),
                 FDR = c(rep(rep(Final_DRs1, times = nM*nTy, each = nY), 3),
                         rep(Final_DRs2, times = nM*nTy, each = nY)),
                 Year = rep(Years, times = nM*nTy*nFt),
                 Mean = rep(NA, nM*nTy*nFt*nY), 
                 SD = rep(NA, nM*nTy*nFt*nY))

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
  Y_sample <- sims_yield[, c(1, 2), , ]
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
  Rel_yield <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
  
  # calculate relative biomass and year arrays after reserve implementation
  for (ty in 1:nTy) {
    for (fdr in 1:nF) {
      for (y in 1:nY) {
        for (sim in 1:num_sims) {
          Rel_biomass[y, ty, fdr, sim] <- B_sample[y + 1, ty, fdr] / 
            B_sample[1, ty, fdr]
          Rel_yield[y, ty, fdr, sim] <- Y_sample[y + 1, ty, fdr] / 
            Y_sample[1, ty, fdr]
        }
      }
    }
  }
  
  ##### fill in data frame with cumulative values #####
  for (ty in 1:nTy) {
    for (fdr in 1:nF) {
      for (y in 1:nY) {
        
        # biomass index
        index_b <- (s - 1)*nM*nTy*nF1*nY +  (ty - 1)*nF*nY + (fdr - 1)*nY + y
        
        # yield index
        index_y <- index_b + nTy*nF*nY
        
        # initialize biomass and yield mean and SD vectors
        bio_mean <- rep(0, num_sims)
        bio_sd <- rep(0, num_sims)
        
        yield_mean <- rep(0, num_sims)
        yield_sd <- rep(0, num_sims)
        
        # iterature through simulations to average values
        for (sim in 1:num_sims) {
          bio_mean[sim] <- mean(Rel_biomass[1:y, ty, fdr, sim]) 
          bio_sd[sim] <- sd(Rel_biomass[1:y, ty, fdr, sim]) 
          
          yield_mean[sim] <- mean(Rel_yield[1:y, ty, fdr, sim]) 
          yield_sd[sim] <- sd(Rel_yield[1:y, ty, fdr, sim]) 
        }
        
        DF$Mean[index_b] <- median(bio_mean)
        DF$SD[index_b]   <- median(bio_sd)
        DF$Mean[index_y] <- median(yield_mean)
        DF$SD[index_y]   <- median(yield_sd)
        
      }
    }
  }
}

# make FDR and Metric factor variables
DF$FDR <- factor(DF$FDR, levels = Final_DRs2)
DF$Metric <- factor(DF$Metric, levels = metrics)
DF$Type <- factor(DF$Type, levels = types)
DF$Species <- factor(DF$Species, levels = Names)

# remove static biomass and yield for low and high estimates
bio_short <- subset(DF, Metric == 'Biomass' & Year == years_to_plot[1])
bio_long <- subset(DF, Metric == 'Biomass' & Year == years_to_plot[2])

yield_short <- subset(DF, Metric == 'Yield' & Year == years_to_plot[1])
yield_long <- subset(DF, Metric == 'Yield' & Year == years_to_plot[2])

# add dummy parameter for x-axis plotting
dif <- 0.2
cr <- 1
br <- 3
lc <- 5
ca <- 7

CR <- c(cr - 3/2*dif, cr - dif/2, cr + dif/2, cr + 3/2*dif)
BR <- c(br - 3/2*dif, br - dif/2, br + dif/2, br + 3/2*dif)
LC <- c(lc - 3/2*dif, lc - dif/2, lc + dif/2, lc + 3/2*dif)
CA <- c(ca - 5/2*dif, ca - 3/2*dif, ca - dif/2, ca + dif/2, 
        ca + 3/2*dif, ca + 5/2*dif)

X <- c(CR, CR, BR, BR, LC, LC, CA, CA)

bio_short$X   <- X
bio_long$X    <- X
yield_short$X <- X
yeild_long$X  <- X

# plotting parameters
og_colors <- rev(viridis(max(c(nF1, nF2)) + 1))
new_colors <- og_colors[2:(nF2 + 1)]

# legend position coordinates
x.legend <- 1.25
y.legend <- 1

##### mean figure #####

# plot cumulative mean biomass year 10
thing1 <- ggplot(bio_short, aes(x = X, y = Mean, color = FDR, 
                                shape = Type)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(stroke = 1, size = 3) +
  ggtitle('(a) Biomass: 10 years') +
  scale_color_manual(values = new_colors) +
  scale_shape_manual(values = c(1, 4)) +
  ylab('Relative Biomass: Mean') +
  scale_x_continuous(breaks = c(cr, br, lc, ca),
                     labels = Names) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = 'none')

# plot cumulative mean yield year 10
thing2 <- ggplot(yield_short, aes(x = X, y = Mean, color = FDR, 
                                  shape = Type)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(stroke = 1, size = 3) +
  ggtitle('(c) Yield: 10 years') +
  scale_color_manual(values = new_colors) +
  scale_shape_manual(values = c(1, 4)) +
  ylab('Relative Yield: Mean') +
  scale_x_continuous(breaks = c(cr, br, lc, ca),
                     labels = Names) +  theme_bw() +
  theme_bw() + 
  theme(axis.title.x = element_blank()) +
  theme(legend.position = 'none')

# plot cumulative mean biomass year 20
thing3 <- ggplot(bio_long, aes(x = X, y = Mean, color = FDR, 
                               shape = Type)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(size = 3, stroke = 1) +
  ggtitle('(b) Biomass: 20 years') +
  scale_color_manual(values = new_colors) +
  scale_shape_manual(values = c(1, 4)) +
  scale_x_continuous(breaks = c(cr, br, lc, ca),
                     labels = Names) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = 'none')

# plot cumulative mean yield year 20
thing4 <- ggplot(yield_long, aes(x = X, y = Mean, color = FDR, 
                                 shape = Type)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point(size = 3, stroke = 1) +
  ggtitle('(d) Yield: 20 years') +
  scale_color_manual(values = new_colors) +
  scale_shape_manual(values = c(1, 4)) +
  scale_x_continuous(breaks = c(cr, br, lc, ca),
                     labels = Names) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(shape = 'Type of \n Control Rule', color = expression('D'[final])) +
  theme(plot.margin = unit(c(0, 100, 0, 0), 'pt')) +
  theme(legend.position = c(x.legend, y.legend))

mean_plot <- (thing1 + thing3) / (thing2 + thing4)

# save results to figures folder
ggsave(mean_plot, filename = paste('fig4_cumulative_mean.png', sep = ''),
       path = paste('C:/Users/Vic/Box Sync/Quennessen_Thesis/MS Thesis/', 
                    'publication manuscript/figures', sep = ''),
       width = png_width, height = png_height)

##### standard deviation figure #####

# plot cumulative SD biomass year 10
thing5 <- ggplot(bio_short, aes(x = X, y = SD, color = FDR, 
                                shape = Type)) +
  geom_point(stroke = 1, size = 3) +
  ggtitle('(a) Biomass: 10 years') +
  scale_color_manual(values = new_colors) +
  scale_shape_manual(values = c(1, 4)) +
  ylab('Relative Biomass: Standard Deviation') +
  scale_x_continuous(breaks = c(cr, br, lc, ca),
                     labels = Names) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = 'none')

# plot cumulative SD yield year 10
thing6 <- ggplot(yield_short, aes(x = X, y = SD, color = FDR, 
                                  shape = Type)) +
  geom_point(stroke = 1, size = 3) +
  ggtitle('(c) Yield: 10 years') +
  scale_color_manual(values = new_colors) +
  scale_shape_manual(values = c(1, 4)) +
  ylab('Relative Yield: Standard Deviation') +
  scale_x_continuous(breaks = c(cr, br, lc, ca),
                     labels = Names) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = 'none')

# plot cumulative SD biomass year 20
thing7 <- ggplot(bio_long, aes(x = X, y = SD, color = FDR, 
                               shape = Type)) +
  geom_point(size = 3, stroke = 1) +
  ggtitle('(b) Biomass: 20 years') +
  scale_color_manual(values = new_colors) +
  scale_shape_manual(values = c(1, 4)) +
  scale_x_continuous(breaks = c(cr, br, lc, ca),
                     labels = Names) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = 'none')

# plot cumulative SD yield year 20
thing8 <- ggplot(yield_long, aes(x = X, y = SD, color = FDR, 
                                 shape = Type)) +
  geom_point(size = 3, stroke = 1) +
  ggtitle('(d) Yield: 20 years') +
  scale_color_manual(values = new_colors) +
  scale_shape_manual(values = c(1, 4)) +
  scale_x_continuous(breaks = c(cr, br, lc, ca),
                     labels = Names) +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(shape = 'Type of \n Control Rule', color = expression('D'[final])) +
  theme(plot.margin = unit(c(0, 100, 0, 0), 'pt')) +
  theme(legend.position = c(x.legend, y.legend))

sd_plot <- (thing5 + thing7) / (thing6 + thing8)

# save results to figures folder
ggsave(sd_plot, filename = paste('fig5_cumulative_SD.png', sep = ''),
       path = paste('C:/Users/Vic/Box Sync/Quennessen_Thesis/MS Thesis/', 
                    'publication manuscript/figures', sep = ''),
       width = png_width, height = png_height)

