library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('recruitment_only.R')

species_list = c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')

mclapply(species_list, recruitment_only, mc.cores = 12, num_sims = 10)
