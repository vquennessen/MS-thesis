source('./cluster.R')

final_DRs <- seq(0.2, 1, by = 0.2)
sims <- 10000

for (i in 1:length(final_DRs)) {
  
  cluster(Species = 'BR_OR_2015', Final_DR = final_DRs[i], num_sims = sims)
  
}
