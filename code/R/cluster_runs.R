source('./base_model.R')
source('./plot_stuff.R')

species <- 'black rockfish 2003'
A <- 5
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
plotting <- F
error <- 0.05

# set numbers of simulations
num_sims <- 1e4

# initialize yield and biomass arrays
sims_yield <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))
sims_biomass <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))
sims_SSB <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))

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
  
}

filepath1 = "../../data/1e4_NM_sims_yield.Rda"
filepath2 = "../../data/1e4_NM_sims_biomass.Rda"
filepath3 = "../../data/1e4_NM_sims_SSB.Rda"

save(NM_sims_yield_1e4, file = filepath1)
save(NM_sims_biomass_1e4, file = filepath2)
save(NM_sims_SSB_1e4, file = filepath3)