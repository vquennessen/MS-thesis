# deterministic runs for overfished scenario (historical FM = 0.23)

# set working directory
setwd('C:/Users/Vic/Documents/Projects/MS-thesis/code')

# load libraries
library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
library(parallel)
source('run_base_model.R')

# species to run
species_list_1 = c('BR_OR_2015')

# final density ratios
Final_DRs_1 <- c(0.6, 0.7, 0.8, 0.9)

# stochastic scenario
Scenario = 'None'

# number of runs
NS1 = 1

# run simulations
lapply(species_list_1, run_base_model, num_sims = NS1, Scenario, Final_DRs_1)
