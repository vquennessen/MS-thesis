source('./code/R/base_model.R')

species <- 'black rockfish 2003'
A <- 5
time1 <- 50
time2 <- 50
allocation <- 'IFD'
R0 <- 1e+5
stochasticity <- T
surveys <- T
fishery_management <- T
fishing <- T
adult_movement <- T

base_model(species, A, time1, time2, allocation, R0, stochasticity, surveys, 
           fishery_management, fishing, adult_movement)