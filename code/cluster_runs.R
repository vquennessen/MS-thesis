source('./cluster.R')

final_DRs <- seq(0.2, 1, by = 0.2)
sims <- 3

for (i in 1:length(final_DRs)) {
  
  cluster(species = 'BR2003', final_DR= final_DRs[i], num_sims = sims)
  
}
