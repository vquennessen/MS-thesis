setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")

source('./base_model.R')
source('./plot_stuff.R')

species <- 'black rockfish 2003'
A <- 5
time1 <- 50
time2 <- 20
CR <- 5
allocation <- 'IFD'
R0 <- 1e+5
stochasticity <- T
surveys <- T
transects <- 24
fishery_management <- T
fishing <- T
adult_movement <- T
plotting <- F

# set numbers of simulations
num_sims <- 1e+05

# initialize yield and biomass arrays
sims_yield <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))
sims_biomass <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))

# run the model for each simulation
for (i in 1:num_sims) {
  
  output <- base_model(species, A, time1, time2, CR, allocation, R0, 
             stochasticity, surveys, transects, fishery_management, 
             fishing, adult_movement, plotting)
  
  # save the relative yield and biomasses for all areas, times after reserve
  # implementation, and control rules
  sims_yield[, , , i] <- output[[1]]
  sims_biomass[, , , i] <- output[[2]]
  
}

save(sims_yield, file = "../../data/1e5_sims_yield.Rda")
save(sims_biomass, file = "../../data/1e5_sims_biomass.Rda")

save(as.data.frame.array(sims_field), file = "../../data/testing_yield.")