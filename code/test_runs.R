library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

species_list = c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')

Scenario = 'Test'

Final_DRs <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 2,
         Scenario, Final_DRs)
