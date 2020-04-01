library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

# species_list = c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')

sims = 2
# mclapply(species_list, run_base_model, mc.cores = 12, num_sims = sims)
# lapply(species_list, run_base_model, num_sims = sims)
lapply('CAB_OR_2019', run_base_model, num_sims = sims)