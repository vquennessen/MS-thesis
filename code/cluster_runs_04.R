remotes::install_github("vquennessen/densityratio")
library(densityratio)
source('./cluster.R')

Species = 'BR_OR_2015'
Final_DR <- 0.4
num_sims <- 1000

cluster(Species, Final_DR, num_sims)
