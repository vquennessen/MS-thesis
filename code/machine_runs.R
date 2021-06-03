setwd('C:/Users/Vic/Documents/Projects/MS-thesis/code')
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

species_list_1 = c('CR_OR_2015', 'LING_OW_2017', 'BR_OR_2015')
species_list_2 = c('CAB_OR_2019')

Final_DRs_1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs_2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

Scenario = 'None'

lapply(species_list_1, run_base_model, num_sims = 1, Scenario, Final_DRs_1)
lapply(species_list_2, run_base_model, num_sims = 1, Scenario, Final_DRs_2)

# Scenario = 'Recruitment'
# 
# lapply(species_list_1, run_base_model, num_sims = 1, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = 1, Scenario, Final_DRs_2)
# 
# Scenario = 'Sampling'
# 
# lapply(species_list_1, run_base_model, num_sims = 1, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = 1, Scenario, Final_DRs_2)
# 
# Scenario = 'Both'
# 
# lapply(species_list_1, run_base_model, num_sims = 1, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = 1, Scenario, Final_DRs_2)
