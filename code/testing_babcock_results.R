rm(list = ls())
setwd('~/Projects/MS-thesis/code')
source('./plot_stuff.R')
devtools::install_git('https://github.com/vquennessen/densityratio.git')
library(densityratio)

# set numbers of simulations
num_sims <- 3

# set arguments
Species = 'BR_CA_2003'
R0 = 1e+5
A = 5
MPAs = c(3)
Time1 = 50
Time2 = 50
Recruitment_mode = 'closed'
Error = 0
Stochasticity = TRUE
Surveys = TRUE
Fishery_management = TRUE
Fishing = TRUE
Transects = 24
Adult_movement = TRUE
Plotting = FALSE
Final_DR = 0.6
Years_sampled = 1
Areas_sampled = 'all'
Ind_sampled = 'all'
Allocation = 'IFD'
Control_rules = c(7:14)
CR <- length(Control_rules)
Plot_individual_runs = F

# total time
TimeT <- Time1 + Time2

# initialize yield and biomass arrays
sims_yield <- array(rep(0, A*TimeT*CR*num_sims), c(A, TimeT, CR, num_sims))
sims_SSB <- array(rep(0, A*TimeT*CR*num_sims), c(A, TimeT, CR, num_sims))

# run the model for each simulation
for (i in 1:num_sims) {
  
  output <- densityratio::base_model(Species, R0, A, MPAs, Time1, Time2, 
                                     Recruitment_mode, Error, Stochasticity, 
                                     Surveys, Fishery_management, Fishing, 
                                     Transects, Adult_movement, Plotting, 
                                     Final_DR, Years_sampled, Areas_sampled, 
                                     Ind_sampled, Allocation, Control_rules)
  
  # save the relative yield and biomasses for all areas, times after reserve
  # implementation, and control rules
  sims_yield[, , , i] <- output[[1]]
  sims_SSB[, , , i] <- output[[3]]

  print(i)
  
}

q <- ifelse(num_sims < 101, num_sims,  paste("1e", log10(num_sims), sep = ''))

filepath1 = paste('../data/', Species, '/', q, "_", Final_DR, "_yield.Rda", 
                  sep = '')
filepath3 = paste('../data/', Species, '/', q, "_", Final_DR, "_SSB.Rda", 
                  sep = '')

save(sims_yield, file = filepath1)
save(sims_SSB, file = filepath3)


