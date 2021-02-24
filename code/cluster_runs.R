library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

# species_list_1 = c('BR_OR_2015', 'LING_OW_2017', 'CR_OR_2015')
# species_list_2 = c('CAB_OR_2019')
# Final_DRs_2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

species_list_1 = c('CR_OR_2015')
Final_DRs_1 <- c(0.6, 0.7, 0.8, 0.9)

########### Recruitment ########################################################

Scenario = 'Recruitment'

mclapply(species_list_1, run_base_model, mc.cores = 12, num_sims = 100,
         Scenario, Final_DRs_1)

########### Sampling ###########################################################

Scenario = 'Sampling'

mclapply(species_list_1, run_base_model, mc.cores = 12, num_sims = 100,
         Scenario, Final_DRs_1)

######### Both #################################################################

Scenario = 'Both'

mclapply(species_list_1, run_base_model, mc.cores = 12, num_sims = 100,
         Scenario, Final_DRs_1)