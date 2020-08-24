setwd('C:/Users/Vic/Documents/Projects/MS-thesis/code')
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

Scenario = 'New.Sampling'
# species_list_1 = c('CR_OR_2015', 'LING_OW_2017', 'BR_OR_2015')
species_list_1 = c('BR_OR_2015')

species_list_2 = c('CAB_OR_2019')

Final_DRs_1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs_2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

lapply(species_list_2, run_base_model, num_sims = 100, Scenario, Final_DRs_2)
lapply(species_list_2, run_base_model, num_sims = 150, Scenario, Final_DRs_2)
lapply(species_list_2, run_base_model, num_sims = 200, Scenario, Final_DRs_2)
lapply(species_list_2, run_base_model, num_sims = 250, Scenario, Final_DRs_2)
lapply(species_list_2, run_base_model, num_sims = 300, Scenario, Final_DRs_2)
lapply(species_list_2, run_base_model, num_sims = 350, Scenario, Final_DRs_2)
lapply(species_list_2, run_base_model, num_sims = 400, Scenario, Final_DRs_2)
lapply(species_list_2, run_base_model, num_sims = 450, Scenario, Final_DRs_2)
lapply(species_list_2, run_base_model, num_sims = 500, Scenario, Final_DRs_2)
