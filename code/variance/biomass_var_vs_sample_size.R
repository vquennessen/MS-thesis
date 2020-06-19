library(ggplot2)

setwd("C:/Users/Vic/Documents/Projects/MS-thesis/code/variance")

num_sims <- 20000
num_medians <- 100
num_variances <- 100

load(paste('../../data/Variance/CAB_OR_2019/', num_sims, '_biomass.Rda', sep = ''))
load(paste('../../data/Variance/CAB_OR_2019/', num_sims, '_yield.Rda', sep = ''))

# set different sample sizes
sample_size <- seq(1e3, num_sims, by = 1e2)

# initialize variance dataframe 
variance_df <- data.frame(Sample.Size = sample_size,
                          B.Variance = rep(NA, length(sample_size)), 
                          Y.Variance = rep(NA, length(sample_size)))

  B_variance <- rep(0, num_variances)
  Y_variance <- rep(0, num_variances)
  
  B_medians1 <- rep(0, num_medians)
  Y_medians1 <- rep(0, num_medians)

for (i in 1:length(sample_size)) {
  
  for (j in 1:num_variances) {
    
    for (k in 1:num_medians) {
      
      indices <- sample(1:num_sims, sample_size[i])
      sampled_biomass <- colSums(sims_biomass[, 21, 1, 1, indices])
      sampled_yield <- sims_yield[21, 1, 1, indices]
      
      # calculate medians 
      B_medians1[k] <- median(sampled_biomass)
      Y_medians1[k] <- median(sampled_yield)
    }
    
    B_variance[j] <- var(B_medians1)
    Y_variance[j] <- var(Y_medians1)
    
  }
  
  variance_df$B.Variance[i] <- median(B_variance)
  variance_df$Y.Variance[i] <- median(Y_variance)
  
}

B <- ggplot(data = variance_df, aes(x = Sample.Size, y = B.Variance)) +
  geom_point() +
  ggtitle('Biomass Variance')

Y <- ggplot(data = variance_df, aes(x = Sample.Size, y = Y.Variance)) +
  geom_line() +
  ggtitle('Yield Variance')

ggsave(filename = paste('Biomass_variance.png', sep = ''), plot = B, 
       path = 'C:/Users/Vic/Box/Quennessen_Thesis/figures/Variance')

ggsave(filename = paste('Yield_variance.png', sep = ''), plot = Y, 
       path = 'C:/Users/Vic/Box/Quennessen_Thesis/figures/Variance')
