library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

species_list = c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')

sims = 40
mclapply(species_list, run_base_model, mc.cores = 12, num_sims = sims)

sims1 = 30
sims2 = 29

mclapply('CAB_OR_2019', run_base_model, mc.cores = 12, num_sims = sims1)
lapply('CAB_OR_2019', run_base_model, num_sims = sims2)

mclapply('LING_OW_2017', run_base_model, mc.cores = 12, num_sims = sims1)
lapply('LING_OW_2017', run_base_model, num_sims = sims2)

mclapply('BR_OR_2015', run_base_model, mc.cores = 12, num_sims = sims1)
lapply('BR_OR_2015', run_base_model, num_sims = sims2)

mclapply('CR_OR_2015', run_base_model, mc.cores = 12, num_sims = sims1)
lapply('CR_OR_2015', run_base_model, num_sims = sims2)
