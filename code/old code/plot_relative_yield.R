# plot relative total yield for static and transient DRs with M

# load any necessary libraries
library(ggplot2)
library(patchwork)
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
num_sims <- 6193
data_folder <- 'Both'
figures_folder <- 'Both/relative'
cluster <- TRUE
png_width <- 7
png_height <- 6
y1 = 0.15
y2 = 1.4
###############################################################################

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')

# set variables
A = 5
MPA = 3
Time2 = 20
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules <- 1:6
types <- c('Static', 'Transient')

# dimensions
PD = 0.25
Error = 0.05
estimates <- c('True', 'Low', 'High')
sample_size = num_sims

nC <- length(Control_rules)
nT <- Time2 + 1
nE <- length(estimates)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1

yield_df <- data.frame(Species = c(rep(Names[1], each = 2*nE*nF1*nT),
                                   rep(Names[2], each = 2*nE*nF1*nT),
                                   rep(Names[3], each = 2*nE*nF1*nT),
                                   rep(Names[4], each = 2*nE*nF2*nT)), 
                       Estimate = c(rep(estimates, each = 2*nF1*nT), 
                                    rep(estimates, each = 2*nF1*nT), 
                                    rep(estimates, each = 2*nF1*nT), 
                                    rep(estimates, each = 2*nF2*nT)),
                       Type = c(rep(types, each = nF1*nT, times = nE), 
                                rep(types, each = nF1*nT, times = nE), 
                                rep(types, each = nF1*nT, times = nE), 
                                rep(types, each = nF2*nT, times = nE)), 
                       FDR = c(rep(Final_DRs1, each = nT, times = 2*nE), 
                               rep(Final_DRs1, each = nT, times = 2*nE), 
                               rep(Final_DRs1, each = nT, times = 2*nE), 
                               rep(Final_DRs2, each = nT, times = 2*nE)), 
                       Time = c(rep(0:Time2, times = 2*nF1*nE), 
                                rep(0:Time2, times = 2*nF1*nE), 
                                rep(0:Time2, times = 2*nF1*nE), 
                                rep(0:Time2, times = 2*nF2*nE)), 
                       Value = rep(0, 2*nE*nFall*nT), 
                       Lower = rep(0, 2*nE*nFall*nT), 
                       Upper = rep(0, 2*nE*nFall*nT))

for (s in 1:length(species_list)) {
  
  # load objects
  if (cluster == TRUE) {
    load(paste('~/Documents/MS-thesis/data/', data_folder, '/', 
               species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
  } else {
    load(paste('~/Projects/MS-thesis/data/', data_folder, '/', 
               species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
  }
  
  Rec_age <- parameters(species_list[s])[[3]]
  Max_age <- parameters(species_list[s])[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  # set nF value for species 
  nF <- ifelse(s == 4, nF2, nF1)  
  
  ##### relative biomass and median, upper, and lower limits  #####
  
  # pull out sample sims
  Y_sample   <- sims_yield[, , , indices]
  
  # initialize relative arrays
  Rel_yield <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
  
  # calculate relative arrays after reserve implementation
  for (cr in 1:nC) {
    for (fdr in 1:nF) {
      for (sim in 1:num_sims) {
        Rel_yield[, cr, fdr, sim] <- Y_sample[, cr, fdr, sim] / 
          Y_sample[1, cr, fdr, sim]
      }
    }
  }
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- Y_lower <- Y_upper <- array(rep(NA, nT*nC*nF), c(nT, nC, nF))
  
  # calculate medians, upper limits, and lower limits  
  for (t in 1:nT) {
    for (cr in 1:nC) {
      for (fdr in 1:nF) {        
        Y_medians[t, cr, fdr] <- median(Rel_yield[t, cr, fdr, ])
        Y_lower[t, cr, fdr] <- quantile(Rel_yield[t, cr, fdr, ], 0.5 - PD)
        Y_upper[t, cr, fdr] <- quantile(Rel_yield[t, cr, fdr, ], 0.5 + PD)
      }
    }
  } 
  
  ##### initialize and fill in data frame #####
  
  for (e in 1:nE) {
    for (type in 1:2) {
      for (fdr in 1:nF) {
        for (t in 1:nT) {
          
          cr <- 2*e - type %% 2
          
          index <- (s - 1)*nE*2*nF1*nT + (e - 1)*2*nF*nT + (type - 1)*nF*nT + (fdr - 1)*nT + t

          yield_df$Value[index] <- Y_medians[t, cr, fdr]
          yield_df$Lower[index] <- Y_lower[t, cr, fdr]
          yield_df$Upper[index] <- Y_upper[t, cr, fdr]
        }
      }
    }
  }
  
}

dfA <- subset(yield_df, Species == Names[1])
dfB <- subset(yield_df, Species == Names[2])
dfC <- subset(yield_df, Species == Names[3])
dfD <- subset(yield_df, Species == Names[4])

# plotting parameters
jitter_height <- 0.015
size1 <- 1.5
size2 <- 0.75
ltype1 <- 2
ltype2 <- 1
colors1 <- c("#00BA38", "#00BFC4", "#619CFF", "#F564E3")
colors2 <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")

# plot it

A <- ggplot(data = dfA, aes(x = Time, y = Value, 
                            color = as.factor(FDR), 
                            linetype = Estimate, 
                            size = Type)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(colors1, 0.05)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_color_manual(values = colors1) +
  scale_size_manual(values = c(size1, size2), guide = 'none') +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') + 
  ggtitle(Names[1]) +
  ylab('Relative Yield') +
  ylim(y1, y2) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

B <- ggplot(data = dfB, aes(x = Time, y = Value, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = as.factor(Type))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(colors1, 0.05)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_color_manual(values = colors1) +
  scale_size_manual(values = c(size1, size2)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ggtitle(Names[2]) +
  ylim(y1, y2) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = 'none')

C <- ggplot(data = dfC, aes(x = Time, y = Value, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = as.factor(Type))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(colors1, 0.05)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_color_manual(values = colors1) +
  scale_size_manual(values = c(size1, size2), guide = 'none') +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') + 
  ggtitle(Names[3]) +
  theme(axis.title.x = element_blank()) +
  ylab('Relative Yield') +
  ylim(y1, y2) +
  theme(legend.position = 'none')

D <- ggplot(data = dfD, aes(x = Time, y = Value, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate), 
                            size = as.factor(Type))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(colors2, 0.05)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height)) +
  scale_color_manual(values = colors2) +
  scale_size_manual(values = c(size1, size2)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
  ggtitle(Names[4]) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y = element_blank()) +
  ylim(y1, y2) +
  xlab('Years since reserve implemented') +
  theme(axis.title.x = element_text(hjust = -4.5)) +
  labs(color = 'FDR', linetype = 'Estimate', size = 'Type') +
  theme(legend.position = c(1.2, 1)) +
  theme(plot.margin = unit(c(0, 80, 0, 0), 'pt')) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         size = guide_legend(order = 3))


##### patch all the figures together #####
patch2 <- (A + B) / (C + D)
thing2 <- patch2 + plot_annotation(title = 'Relative Yield')

if (cluster == TRUE) {
  ggsave(thing2, filename = 'M_relative_yield.png',
         path = paste('~/Documents/MS-thesis/figures/', figures_folder, sep = ''),
         width = png_width, height = png_height)
} else {
  ggsave(thing2, filename = 'M_relative_yield.png',
         path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                      figures_folder, sep = ''),
         width = png_width, height = png_height)
}
