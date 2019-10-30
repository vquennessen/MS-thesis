library(viridis)

setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")

load("../../data/1e4_sims_SSB.Rda")

num_sims <- 1e4

y1 <- 1.5
y2 <- 2.5

# create empty plot
par(mar=c(4.5,5,3,1))
plot(1, type = 'l',                          # make an empty line graph
     main = 'Variance of SSB',             # title of plot
     ylab = 'Sample Variance',               # axis labels
     xlab = 'Sample Size', 
     xaxt = 'n', 
     yaxt = 'n',
     xlim = c(0, num_sims),
     ylim = c(y1, y2), 
     cex.main = 2, 
     cex.lab = 2)

# set specific y-axis
ytick <- seq(y1, y2, by = (y2 - y1)/2)            # set y axis tick marks
axis(side = 2,                               # specify y axis
     at = ytick,                             # apply tick marks
     labels = T,                             # apply appropriate labels
     las = 0,                                # set text horizontal
     cex.axis = 1.5)

# set specific x-axis
xtick <- seq(0, num_sims, by = num_sims/4)   # set x axis tick marks
axis(side = 1,                               # specify x axis
     at = xtick,                             # apply tick marks
     labels = T,                             # apply appropriate labels
     las = 1,                                # set text horizontal    
     cex.axis = 1.5)

num_runs <- 10
time2 <- 20
color <- viridis(num_runs)

# set different sample sizes
sample_size <- seq(1e2, num_sims, by = 1e2)

# initialize variance arrays 
SSB_vars <- rep(NA, length(sample_size))

for (j in 1:num_runs) {
  
  for (i in 1:length(sample_size)) {
    
    indices <- sample(1:num_sims, sample_size[i])
    sampled_SSB <- sims_SSB[, , , indices]
    
    # calculate medians 
    SSB_vars[i] <- var(sampled_SSB[1, time2, 1, ])
    
  }
  
  # plot one run through various sample sizes
  lines(sample_size, SSB_vars, type = 'l', col = color[j])
}
