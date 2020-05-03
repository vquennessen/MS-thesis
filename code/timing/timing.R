#load packages
library(generics)
library(ggplot2)

# read in timing spreadsheet
timing <- read.csv('~/Projects/MS-thesis/code/timing.csv', header = T)

# convert time difference to minutes
timing$minutes <- round(as.difftime(as.character(timing$time_difference), 
                                    format = '%H:%M:%S', units = 'mins'), 2)

# divide by sims
timing$time_per_sim <- timing$minutes / timing$num_sims

# visualize time per sim by function
ggplot(data = timing, aes(x = species, y = time_per_sim, fill = function.)) +
  geom_boxplot()
