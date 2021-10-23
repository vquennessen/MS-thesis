# proportion analysis

# set working directory
# setwd("~/Projects/MS-thesis/code")

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

ggsave(fig2, file = 'fig7.png', 
       path = '~/Documents/MS-thesis/figures',
       width = 5, height = 6)
