library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

Scenario = 'Variance'

species_list = c('CAB_OR_2019')

Final_DRs <- c(0.9)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 1050,
         Scenario, Final_DRs)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 2050,
         Scenario, Final_DRs)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 3050,
         Scenario, Final_DRs)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 4050,
         Scenario, Final_DRs)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 950,
         Scenario, Final_DRs)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 1950,
         Scenario, Final_DRs)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 2950,
         Scenario, Final_DRs)

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 3950,
         Scenario, Final_DRs)
