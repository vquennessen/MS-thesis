library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('no_stochasticity.R')

species_list = c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')

for (s in 1:length(species_list)) {
  no_stochasticity(Species = species_list[s], num_sims = 2)
}
