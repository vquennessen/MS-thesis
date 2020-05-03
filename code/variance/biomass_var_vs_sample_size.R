library(ggplot2)

setwd("C:/Users/Vic/Documents/Projects/MS-thesis/code/variance")

num_sims <- 21425
num_runs <- 10
times <- c(5, 10, 15, 20)

load(paste('../../data/variance_runs/', num_sims, '_biomass.Rda', sep = ''))
load(paste('../../data/variance_runs/', num_sims, '_yield.Rda', sep = ''))

# set different sample sizes
sample_size <- seq(1e2, num_sims, by = 1e2)


for (t in 1:length(times)) {
  
  # initialize variance dataframe 
  variance_df <- data.frame(Run = rep(1:num_runs, each = length(sample_size)),
                            Sample.Size = rep(sample_size, times = num_runs),
                            B.Variance = rep(NA, length(sample_size)*num_runs), 
                            Y.Variance = rep(NA, length(sample_size)*num_runs))
  
  for (j in 1:num_runs) {
    
    for (i in 1:length(sample_size)) {
      
      indices <- sample(1:num_sims, sample_size[i])
      sampled_biomass <- sims_biomass[, , , , indices]
      sampled_yield <- sims_yield[, , , , indices]
      
      # calculate variance 
      index <- (j - 1)*length(sample_size) + i
      variance_df$B.Variance[index] <- var(sampled_biomass[1, times[t], 1, ])
      variance_df$Y.Variance[index] <- var(sampled_yield[1, times[t], 1, ])
    }
    
  }
  
  B <- ggplot(data = variance_df, aes(x = Sample.Size, y = B.Variance, 
                                      color = as.factor(Run))) +
    geom_line() +
    ggtitle('Biomass Variance') +
    labs(color = 'Run')
  
  Y <- ggplot(data = variance_df, aes(x = Sample.Size, y = Y.Variance, 
                                      color = as.factor(Run))) +
    geom_line() +
    ggtitle('Yield Variance') +
    labs(color = 'Run')
  
  ggsave(paste('Biomass_variance_year', times[t], '.png', sep = ''), B, 
         path = 'C:/Users/Vic/Google Drive/OSU/Thesis/figures/variance')
  
  ggsave(paste('Yield_variance_year', times[t], '.png', sep = ''), Y, 
         path = 'C:/Users/Vic/Google Drive/OSU/Thesis/figures/variance')
  
}