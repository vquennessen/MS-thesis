setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")

num_sims <- 1e4
species <- 'BR2003'
final_DR <- 0.8
factor <- 'yield'

filepath <- paste('../../data/1e', log10(num_sims), '_', species, '_', final_DR, 
                  '_', factor, '.Rda', sep = '')

load(filepath)

y1 <- 0
y2 <- 5

# create empty plot
par(mar=c(4.5,5,3,1))
plot(1, type = 'l',                          # make an empty line graph
     main = 'Variance vs. Sample Size',      # title of plot
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

# set different sample sizes
sample_size <- seq(1e2, num_sims, 1e2)

# initialize variance arrays 
vars <- array(rep(NA, length(sample_size*num_runs)), c(sample_size, num_runs))

for (j in 1:num_runs) {
  
  for (i in 1:length(sample_size)) {
    
    indices <- sample(1:num_sims, sample_size[i])
    sampled <- sims_yield[, , , indices]
    
    # calculate variance of yield at time2 across all simulations 
    vars[i] <- var(sampled[1, time2, 1, indices])
    
  }
  
  # plot one run through various sample sizes
  points(sample_size, vars[, j])
  
}