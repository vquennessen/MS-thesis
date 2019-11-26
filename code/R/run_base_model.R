rm(list = ls())

setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")

source('./base_model.R')
source('./plot_stuff.R')

species <- 'black rockfish 2003'
A <- 5
MPAs <- c(3)
time1 <- 50
time2 <- 20
CR <- 6
allocation <- 'IFD'
R0 <- 1e+5
stochasticity <- T
surveys <- T
transects <- 24
fishery_management <- T
fishing <- T
adult_movement <- T
error <- 0.05
final_DR <- 0.8
plotting <- T
plot_individual_runs <- F

# set numbers of simulations
num_sims <- 1e2

# initialize yield and biomass arrays
sims_yield <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))
sims_biomass <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))
sims_SSB <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))
sims_DR <- array(rep(0, time2*CR*num_sims), c(time2, CR, num_sims))

# run the model for each simulation
for (i in 1:num_sims) {
  
  output <- base_model(species, A, time1, time2, CR, allocation, R0, 
             stochasticity, surveys, transects, fishery_management, 
             fishing, adult_movement, plotting, error)
  
  # save the relative yield and biomasses for all areas, times after reserve
  # implementation, and control rules
  sims_yield[, , , i] <- output[[1]]
  sims_biomass[, , , i] <- output[[2]]
  sims_SSB[, , , i] <- output[[3]]
  sims_DR[, , i] <- output[[4]]
  
}

filepath1 = "../../data/1e2_NM_sims_yield.Rda"
filepath2 = "../../data/1e2_NM_sims_biomass.Rda"
filepath3 = "../../data/1e2_NM_sims_SSB.Rda"
filepath4 = "../../data/1e2_NM_sims_DR.Rda"

save(sims_yield, file = filepath1)
save(sims_biomass, file = filepath2)
save(sims_SSB, file = filepath3)
save(sims_DR, file = filepath4)

plot_stuff(filepath1, filepath2, filepath3, filepath4, A, time2, CR, num_sims, 
           sample_size = num_sims, PD = 0.25, plot_individual_runs)
