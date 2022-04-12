# plot difference in total yield for static vs. transient DRs

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
figures_folder <- 'Both/difference'
cluster <- TRUE
png_width <- 6
png_height <- 5
y1 <- -0.75
y2 <- 1
y1.1 <- -0.25
y2.1 <- 2.25
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

# dimensions
sample_size = num_sims
PD = 0.25
plot_individual_runs = FALSE
Error = 0.05
estimates <- c('True', 'Low', 'High')
ENM = 2

nT <- Time2 + 1
nE <- length(estimates)
nC <- length(Control_rules)
nS <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFall <- nF2 + 3*nF1

yield_df <- data.frame(Species = c(rep(Names[1], each = nE*nF1*nT),
                                   rep(Names[2], each = nE*nF1*nT),
                                   rep(Names[3], each = nE*nF1*nT),
                                   rep(Names[4], each = nE*nF2*nT)), 
                       Estimate = c(rep(estimates, each = nF1*nT), 
                                    rep(estimates, each = nF1*nT), 
                                    rep(estimates, each = nF1*nT), 
                                    rep(estimates, each = nF2*nT)), 
                       FDR = c(rep(Final_DRs1, each = nT, times = nE), 
                               rep(Final_DRs1, each = nT, times = nE), 
                               rep(Final_DRs1, each = nT, times = nE), 
                               rep(Final_DRs2, each = nT, times = nE)), 
                       Time = c(rep(0:Time2, times = nF1*nE), 
                                rep(0:Time2, times = nF1*nE), 
                                rep(0:Time2, times = nF1*nE), 
                                rep(0:Time2, times = nF2*nE)), 
                       Difference = rep(0, nE*nFall*nT), 
                       Lower = rep(0, nE*nFall*nT), 
                       Upper = rep(0, nE*nFall*nT))

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
  M <- parameters(species_list[s])[[2]]
  Nat_mortality <- c(M - Error, M, M + Error)
  NM <- length(Nat_mortality)
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  # set nF value for species 
  nF <- ifelse(s == 4, nF2, nF1)  
 
  # pull out sample sims
  Y_sample   <- sims_yield[, , , indices]
  
  # initialize relative arrays
  Rel_yield <- array(rep(0, nT*nC*nF*num_sims), c(nT, nC, nF, num_sims))
  
  # calculate relative arrays after reserve implementation
  for (cr in 1:nC) {
    for (fdr in 1:nF) {
      for (sim in indices) {
        Rel_yield[, cr, fdr, sim] <- Y_sample[, cr, fdr, sim] / Y_sample[1, cr, fdr, sim]
      }
    }
  }
  ##### initialize difference array #####
  Difference <- array(rep(0, nT*nE*nF*num_sims), c(nT, nE, nF, num_sims))
  
  ##### calculate differences between transient and static DRs #####
    for (t in 1:nT) {
      for (e in 1:nE) {
        for (fdr in 1:nF) {     
          for (sim in 1:num_sims) {
            diff <- Rel_yield[t, 2*e, fdr, sim] - Rel_yield[t, 2*e - 1, fdr, sim]
            Difference[t, e, fdr, sim] <- diff / Rel_yield[t, 2*e - 1, fdr, sim]
          }
        }
      }
    }
  
  ##### fill in data frames with median and quantile values #####
    for (e in 1:nE) {
      for (fdr in 1:nF) {
        for (t in 1:nT) {
        index <- (s - 1)*nE*nF1*nT + (e - 1)*nF*nT + (fdr - 1)*nT + t
        yield_df$Difference[index] <- median(Difference[t, e, fdr, ])
        yield_df$Lower[index] <- quantile(Difference[t, e, fdr, ], 0.5 - PD)
        yield_df$Upper[index] <- quantile(Difference[t, e, fdr, ], 0.5 + PD)

      }
    }
  }

}

dfA <- subset(yield_df, Species == Names[1])
dfB <- subset(yield_df, Species == Names[2])
dfC <- subset(yield_df, Species == Names[3])
dfD <- subset(yield_df, Species == Names[4])

# plotting parameters
jitter_height <- 0.005
if (s != 4) {
  colors <- c("#00BA38", "#00BFC4", "#619CFF", "#F564E3")
} else {
  colors <- c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3")
}

# plot it
A <- ggplot(data = dfA, aes(x = Time, y = Difference, 
                           color = as.factor(FDR), 
                           linetype = as.factor(Estimate))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(c(colors), 0.15)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height2)) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  labs(color = 'FDR', linetype = 'Estimate of M') +
  ggtitle(Names[1]) +
  ylab('Difference') +
  ylim(y1, y2) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank())

B <- ggplot(data = dfB, aes(x = Time, y = Difference, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(c(colors), 0.15)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height2)) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  labs(color = 'FDR', linetype = 'Estimate of M') +
  ggtitle(Names[2]) +
  ylim(y1, y2) +
  theme(axis.text = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = 'none')

C <- ggplot(data = dfC, aes(x = Time, y = Difference, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(c(colors), 0.15)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height2)) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  labs(color = 'FDR', linetype = 'Estimate of M') +
  ggtitle(Names[3]) +
  theme(axis.title.x = element_blank()) +
  ylab('Difference') +
  ylim(y1, y2) +
  theme(legend.position = 'none')

D <- ggplot(data = dfD, aes(x = Time, y = Difference, 
                            color = as.factor(FDR), 
                            linetype = as.factor(Estimate))) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(FDR), 
                  colour = NA), show.legend = FALSE) +
  scale_fill_manual(values = alpha(c(colors), 0.15)) +
  geom_line(position = position_jitter(w = 0, h = jitter_height2)) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black') + 
  ggtitle(Names[4]) +
  theme(axis.title.y = element_blank()) +
  ylim(y1.1, y2.1) +
  xlab('Years since reserve implemented') +
  theme(axis.title.x = element_text(hjust = 5)) +
  theme(legend.position = c(1.25, 1.2)) +
  theme(plot.margin = unit(c(0, 70, 0, 0), 'pt')) +
  labs(color = 'FDR', linetype = 'Estimate \n of M')  +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))

##### patch all the figures together #####
patch2 <- (A + B) / (C + D)
thing2 <- patch2 + plot_annotation(
  title = 'Relative Yield (Transient - Static DRCR)')

if (cluster == TRUE) {
  ggsave(thing2, filename = 'M_yield_difference.png',
         path = paste('~/Documents/MS-thesis/figures/', figures_folder, sep = ''),
         width = png_width, height = png_height)
} else {
  ggsave(thing2, filename = 'M_yield_difference.png',
         path = paste('C:/Users/Vic/Box/Quennessen_Thesis/figures/', 
                      figures_folder, sep = ''),
         width = png_width, height = png_height)
}

