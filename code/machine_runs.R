setwd('C:/Users/Vic/Documents/Projects/MS-thesis/code')
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

species_list_1 = c('CR_OR_2015_SSS', 'LING_OW_2017_SSS', 'BR_OR_2015_SSS')
species_list_2 = c('CAB_OR_2019_SSS')

Final_DRs_1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs_2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

NS1 <- 1
# NS2 <- 25
# NS3 <- 20
# NS4 <- 15

Scenario = 'None'

lapply(species_list_1, run_base_model, num_sims = NS1, Scenario, Final_DRs_1)
lapply(species_list_2, run_base_model, num_sims = NS1, Scenario, Final_DRs_2)

# lapply(species_list_1, run_base_model, num_sims = NS2, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = NS2, Scenario, Final_DRs_2)
# 
# lapply(species_list_1, run_base_model, num_sims = NS3, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = NS3, Scenario, Final_DRs_2)
# 
# lapply(species_list_1, run_base_model, num_sims = NS4, Scenario, Final_DRs_1)
# lapply(species_list_2, run_base_model, num_sims = NS4, Scenario, Final_DRs_2)
