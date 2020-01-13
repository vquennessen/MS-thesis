rm(list = ls())

setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")

source('./base_model.R')
source('./plot_stuff.R')
source('./transient_DR.R')
library(densityratio)

Species <- 'BR2003'
Stochasticity <- T
Surveys <- T
Fishery_management <- T
Fishing <- T
Adult_movement <- T
Final_DR <- 0.2
Plotting <- T
Plot_individual_runs <- F
Years_sampled <- 1

# set numbers of simulations
num_sims <- 3

# initialize yield and biomass arrays
sims_yield <- array(rep(0, A*(time2 + 1)*CR*num_sims), 
                    c(A, (time2 + 1), CR, num_sims))
sims_biomass <- array(rep(0, A*(time2 + 1)*CR*num_sims), 
                      c(A, (time2 + 1), CR, num_sims))
sims_SSB <- array(rep(0, A*(time2 + 1)*CR*num_sims), 
                  c(A, (time2 + 1), CR, num_sims))
sims_DR <- array(rep(0, (time2 + 1)*CR*num_sims), 
                 c((time2 + 1), CR, num_sims))
y_DR <- array(rep(0, (time2 + 1)*num_sims), c((time2 + 1), num_sims))

# run the model for each simulation
for (i in 1:num_sims) {
  
  output <- base_model(Species, Stochasticity, Surveys, Fishery_management, 
                       Fishing, Adult_movement, Plotting, Final_DR, 
                       Years_sampled)
  
  # save the relative yield and biomasses for all areas, times after reserve
  # implementation, and control rules
  sims_yield[, , , i] <- output[[1]]
  sims_biomass[, , , i] <- output[[2]]
  sims_SSB[, , , i] <- output[[3]]
  sims_DR[, , i] <- output[[4]]
  y_DR[, i] <- output[[5]]
  
}

q <- ifelse(num_sims < 101, num_sims,  paste("1e", log10(num_sims), sep = ''))

filepath1 = paste('../../data/', species, '/', q, "_", final_DR, "_yield.Rda", 
                  sep = '')
filepath2 = paste('../../data/', species, '/', q, "_", final_DR, "_biomass.Rda", 
                  sep = '')
filepath3 = paste('../../data/', species, '/', q, "_", final_DR, "_SSB.Rda", 
                  sep = '')
filepath4 = paste('../../data/', species, '/', q, "_", final_DR, "_DR.Rda", 
                  sep = '')

save(sims_yield, file = filepath1)
save(sims_biomass, file = filepath2)
save(sims_SSB, file = filepath3)
save(sims_DR, file = filepath4)

plot_stuff(filepath1, filepath2, filepath3, filepath4, A, Time2, CR, num_sims, 
           sample_size = num_sims, PD = 0.25, Plot_individual_runs, 
           y_DR[, num_sims], Species, Final_DR)