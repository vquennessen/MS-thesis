rm(list = ls())
setwd('~/Projects/MS-thesis/code')
source('./plot_stuff.R')
devtools::install_git('https://github.com/vquennessen/densityratio.git')
library(densityratio)

# set numbers of simulations
num_sims <- 250

# set arguments
Species = 'BR_OR_2015'
R0 = 1e+5
A = 5
MPAs = c(3)
Time1 = 50
Time2 = 20
Recruitment_mode = 'pool'
Error = 0.05
Stochasticity = TRUE
Surveys = TRUE
Fishery_management = TRUE
Fishing = TRUE
Transects = 24
Adult_movement = TRUE
Plotting = FALSE
Final_DR = 0.2
Years_sampled = 1
Areas_sampled = 'all'
Ind_sampled = 'all'
Allocation = 'IFD'
Control_rules = c(1:6)
CR <- length(Control_rules)
Plot_individual_runs <- F
BM <- FALSE

y_DR <- densityratio::transient_DR(Time1 = 50, TimeT = 70, Final_DR,
                                   Nat_mortality = c(0.09, 0.14, 0.19), nm = 2)

# total time
TimeT <- Time1 + Time2

# initialize yield and biomass arrays
sims_yield <- array(rep(0, 2*TimeT*CR*1*num_sims), c(2, TimeT, CR, 1, num_sims))
sims_biomass <- array(rep(0, 3*TimeT*CR*1*num_sims), c(3, TimeT, CR, 1, num_sims))
sims_SSB <- array(rep(0, 3*TimeT*CR*1*num_sims), c(3, TimeT, CR, 1, num_sims))
sims_DR <- array(rep(0, (Time2 + 1)*CR*num_sims), c(Time2 + 1, CR, num_sims))
sims_N <- array(rep(0, 38*3*TimeT*CR*1*num_sims), c(38, 3, TimeT, CR, 1, num_sims))

# run the model for each simulation
for (i in 1:num_sims) {
  
  output <- densityratio::base_model(Species, R0, A, MPAs, Time1, Time2, 
                                     Recruitment_mode, Error, Stochasticity, 
                                     Surveys, Fishery_management, Fishing, 
                                     Transects, Adult_movement, Plotting, 
                                     Final_DR, Years_sampled, Areas_sampled, 
                                     Ind_sampled, Allocation, BM, Control_rules)
  
  # save the relative yield and biomasses for all areas, times after reserve
  # implementation, and control rules
  sims_yield[, , , , i] <- output[[1]]
  sims_biomass[, , , , i] <- output[[2]]
  sims_SSB[, , , , i] <- output[[3]]
  sims_DR[, , i] <- output[[4]]
  sims_N[, , , , , i] <- output[[5]]

  print(i)
  
}

q <- ifelse(num_sims < 1000, num_sims,  paste("1e", log10(num_sims), sep = ''))

filepath1 = paste('home/quennessenv/ExpanDrive/Box/data/', Species, '/', q, "_", 
                  Final_DR, "_yield.Rda", sep = '')
filepath2 = paste('home/quennessenv/ExpanDrive/Box/data/', Species, '/', q, "_", 
                  Final_DR, "_biomass.Rda", sep = '')
filepath3 = paste('home/quennessenv/ExpanDrive/Box/data/', Species, '/', q, "_", 
                  Final_DR, "_SSB.Rda", sep = '')
filepath4 = paste('home/quennessenv/ExpanDrive/Box/data/', Species, '/', q, "_", 
                  Final_DR, "_DR.Rda", sep = '')
filepath5 = paste('home/quennessenv/ExpanDrive/Box/data/', Species, '/', q, '-', 
                  Final_DR, '_N.Rda', sep = '')

save(sims_yield, file = filepath1)
save(sims_biomass, file = filepath2)
save(sims_SSB, file = filepath3)
save(sims_DR, file = filepath4)
save(sims_N, file = filepath5)

plot_stuff(filepath1, filepath2, filepath3, filepath4, filepath5, A, Time2, CR, 
           num_sims, sample_size = num_sims, PD = 0.25, Plot_individual_runs, 
           y_DR, Species, Final_DR)
