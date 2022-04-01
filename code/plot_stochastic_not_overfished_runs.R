# plot stochastic (not overfished) runs for BR OR 2015

##### figure 6/S11: stochastic vs deterministic, transient vs. static ##########

# load any necessary libraries
library(ggplot2)
library(patchwork)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(egg)
library(viridis)

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

# species to compare
species_list <- c('BR_OR_2015')
Names <- c('Black Rockfish')
titles <- c('figS11_', 'fig6_')

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

for (v in 1:2) {
  
  if (v == 1) { 
    FDRs <- list(c(0.6, 0.9))
    
  } else if (v == 2) {
    FDRs <- c(0.7)
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
    
    ##### plotting parameters #####
    jitter_height <- 0
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
      
      final_plot <- ggplot(data = new_DF, aes(x = Year, y = Value, color = FDR, 
                                              linetype = Type)) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = FDR, colour = NA), 
                    show.legend = FALSE) +
        scale_fill_manual(values = alpha(new_colors, alfa)) +
        geom_line(position = position_jitter(w = 0, h = jitter_height)) +
        scale_color_manual(values = new_colors) +
        geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
        theme_bw() +
        theme(axis.title.y = element_blank()) +
        facet_grid(Metric ~ Scenario, scales = 'free', switch = 'y') +
        xlab('Years since reserve implemented') +
        guides(color = guide_legend(order = 1), 
               linetype = guide_legend(order = 2)) +
        labs(color = expression('D'[final]), 
             linetype = 'Type of \n Control Rule')
      
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
      final_plot <- ggplot(data = combined_DF, 
                           aes(x = Year, y = Value, color = Type, 
                               linetype = Source)) +
        geom_hline(yintercept = 1, linetype = 'dashed', color = 'black') +
        geom_ribbon(data = new_DF, aes(ymin = Lower, ymax = Upper, fill = Type,
                                       colour = NA), show.legend = FALSE) +
        scale_fill_manual(values = alpha(new_colors, alfa)) +
        geom_line(position = position_jitter(w = 0, h = jitter_height)) +
        scale_color_manual(values = new_colors) +
        theme_bw() +
        theme(axis.title.y = element_blank()) +
        facet_grid(Metric ~ Scenario, scales = 'free', switch = 'y') +
        xlab('Years since reserve implemented') +
        guides(color = guide_legend(order = 1), 
               linetype = guide_legend(order = 2)) +
        labs(color = 'Type of \n Control Rule', linetype = 'Source')
      
    }
    
    final_plot <- tag_facet(final_plot)
    final_plot <- final_plot +
      theme(strip.text = element_text(), strip.background = element_rect())
    
    # figure titles index
    fti <- (v - 1)*length(species_list) + s
    
    # save plot
    ggsave(final_plot,
           filename = paste(titles[fti], Names[s], '_not_overfished.png', sep = ''),
           path = '~/Documents/MS-thesis/figures/not_overfished',
           width = png_width, height = png_height)
  }
  
}

################################################################################

##### figure 7: proportions across all species #################################

# load libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(egg)

# source function to create proportion_transient_greater.csv
source('extract_proportions_transient_greater.R')

# load data
data <- read.csv("proportion_transient_greater.csv", header = TRUE)

# calculate means by grouping species by FDR value
basic_FDR_year <- data %>%
  group_by(Metric, Species, FDR) %>%
  dplyr::summarize(Mean = mean(Value, na.rm=TRUE))

basic_FDR_year

# metrics
metrics <- c('Biomass', 'Yield', 'Effort')
metric_names <- c('Relative Biomass', 'Relative Yield', 'Relative Effort')
basic_FDR_year$Metric <- factor(basic_FDR_year$Metric,levels = metrics, 
                                labels = metric_names)

# colors
n <- length(unique(basic_FDR_year$Species))
og_colors <- viridis(n + 1)
colors <- og_colors[1:n]

# figure width and height
png_width <- 6

# make plot
fig1 <- ggplot(data = basic_FDR_year, 
               aes(x = FDR, y = Mean, color = Species, shape = Species)) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', size = 0.75) +
  geom_line(size = 0.75, alpha = 0.75) +
  geom_point(size = 3, alpha = 0.75) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c(15, 16, 17, 18)) + 
  facet_grid(rows = vars(Metric), scales = 'free', switch = 'y') +
  ylab('Mean Proportion of Simulations') +
  xlab(expression('D'[final])) +
  theme_bw()

fig2 <- tag_facet(fig1, vjust = 2.5) +
  theme(strip.text = element_text(), strip.background = element_rect())

ggsave(fig2, file = 'fig7_not_overfished.png', 
       path = '~/Documents/MS-thesis/figures/not_overfished',
       width = 5, height = 6)

