library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

species_list = c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')

mclapply(species_list, run_base_model, mc.cores = 12, num_sims = 72)

lapply('CR_OR_2015', run_base_model, num_sims = 52)
lapply('CR_OR_2015', run_base_model, num_sims = 100)
lapply('CR_OR_2015', run_base_model, num_sims = 100)

lapply('CAB_OR_2019', run_base_model, num_sims = 95)
lapply('CAB_OR_2019', run_base_model, num_sims = 96)
lapply('CAB_OR_2019', run_base_model, num_sims = 97)
lapply('CAB_OR_2019', run_base_model, num_sims = 98)
lapply('CAB_OR_2019', run_base_model, num_sims = 99)
lapply('CAB_OR_2019', run_base_model, num_sims = 101)
lapply('CAB_OR_2019', run_base_model, num_sims = 102)
lapply('CAB_OR_2019', run_base_model, num_sims = 103)
lapply('CAB_OR_2019', run_base_model, num_sims = 104)
lapply('CAB_OR_2019', run_base_model, num_sims = 105)
