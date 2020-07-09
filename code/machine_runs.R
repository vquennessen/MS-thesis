library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

Scenario = 'Sampling.HighT'
species_list_1 = c('CR_OR_2015', 'LING_OW_2017', 'BR_OR_2015')
# species_list_1 = c('CR_OR_2015', 'BR_OR_2015')

species_list_2 = c('CAB_OR_2019')

Final_DRs_1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs_2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

lapply(species_list_1, run_base_model, num_sims = 200, Scenario, Final_DRs_1)
lapply(species_list_2, run_base_model, num_sims = 200, Scenario, Final_DRs_2)

# lapply(species_list_1, run_base_model, num_sims = 10, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = 10, Scenario, Final_DRs_2)
# 
# lapply(species_list_1, run_base_model, num_sims = 10, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = 10, Scenario, Final_DRs_2)


# lapply(species_list_1, run_base_model, num_sims = 10, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = 10, Scenario, Final_DRs_2)
# 
# lapply(species_list_1, run_base_model, num_sims = 10, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = 10, Scenario, Final_DRs_2)