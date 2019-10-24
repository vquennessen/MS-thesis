rm(list = ls())

source('./code/R/base_model.R')

species <- 'black rockfish 2015'
A <- 5
time1 <- 50
time2 <- 50
CR <- 8
allocation <- 'IFD'
R0 <- 1e+5
stochasticity <- T
surveys <- T
transects <- 24
fishery_management <- T
fishing <- T
adult_movement <- T

# set numbers of simulations
sims <- c(1, 2, 10, 100)

data <- as.data.frame(sims)
data$time <- rep(0, length(sims))

for (j in 1:length(sims)) {
  
  start <- Sys.time()
  
  num_sims <- sims[j]
  
  # initialize yield and biomass arrays
  sims_yield <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))
  sims_biomass <- array(rep(0, A*time2*CR*num_sims), c(A, time2, CR, num_sims))
  
  # run the model for each simulation
  for (i in 1:num_sims) {
    
    output <- base_model(species, A, time1, time2, CR, allocation, R0, 
                         stochasticity, surveys, transects, fishery_management, 
                         fishing, adult_movement)
    
    # save the relative yield and biomasses for all areas, times after reserve 
    # implementation, and control rules
    sims_yield[, , , i] <- output[[1]]
    sims_biomass[, , , i] <- output[[2]]
    
  }
  
  # write the simulation results to .csv files
  write.table(sims_yield, row.names = F, col.names = F, append = F,
              file = "C:/Users/Vic/Documents/Projects/DensityRatio/data/testing_yield.csv")
  
  write.table(sims_biomass, row.names = F, col.names = F, append = F,
              file = "C:/Users/Vic/Documents/Projects/DensityRatio/data/testing_biomass.csv")
  
  end <- Sys.time()
  
  data$time[j] <- difftime(end, start, units = 'mins')
}

data