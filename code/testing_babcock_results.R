rm(list = ls())
setwd('~/Projects/MS-thesis/code')
devtools::install_git('https://github.com/vquennessen/densityratio.git')
library(densityratio)

# set numbers of simulations
num_sims <- 5

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
Years_sampled = NULL
Areas_sampled = NULL
Ind_sampled = NULL
Allocation = 'IFD'
BM = TRUE
Control_rules = c(1:8)
CR <- length(Control_rules)
Plot_individual_runs = FALSE

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
                                     Ind_sampled, Allocation, BM, Control_rules)
  
  # save the relative yield and biomasses for all areas, times after reserve
  # implementation, and control rules
  sims_yield[, , , i] <- output[[1]]
  sims_SSB[, , , i] <- output[[3]]

  print(i)
  
} 
  
  # initialize arrays for sums of yield and SSB
  Y <- array(rep(0, TimeT*CR*num_sims), c(TimeT, CR, num_sims))
  SSB <- array(rep(0, TimeT*CR*num_sims), c(TimeT, CR, num_sims))
  
  for (t in 1:TimeT) {
    for (cr in 1:CR) {
      for (sim in 1:num_sims) {
        Y[t, cr, sim] <- sum(sims_yield[, t, cr, sim])
        SSB[t, cr, sim] <- sum(sims_SSB[, t, cr, sim])
      }
    }
  }
  
  # starting yield and biomass
  start_Y <- array(rep(0, CR*num_sims), c(CR, num_sims))
  start_SSB <- array(rep(0, CR*num_sims), c(CR, num_sims))
  
  for (cr in 1:CR) {
    for (sim in 1:num_sims) {
      start_Y[cr, sim] <- sum(sims_yield[c(1, 2, 4, 5), Time1, cr, sim])
      start_SSB[cr, sim] <- sum(sims_SSB[c(1, 2, 4, 5), Time1, cr, sim])
    }
  }
  
  # relative yield and SSB after reserve implementation
  rel_Y <- array(rep(0, (Time2 + 1)*CR*num_sims), c(Time2 + 1, CR, num_sims))
  rel_SSB <- array(rep(0, (Time2 + 1)*CR*num_sims), c(Time2 + 1, CR, num_sims))
  
  for (t in Time1:TimeT) {
    for (cr in 1:CR) {
      for (sim in 1:num_sims) {
        rel_Y[t - Time1 + 1, cr, sim] <- Y[t, cr, sim]/start_Y[cr, sim]
        rel_SSB[t - Time1 + 1, cr, sim] <- SSB[t, cr, sim]/start_SSB[cr, sim]
      }
    }
  }
  
  # pull out only every 5 years
  indices <- seq(Time1 + 10, TimeT, by = 10)
  ind <- length(indices)
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_means <- array(rep(NA, ind*CR), c(ind, CR))
  SSB_means <- array(rep(NA, ind*CR), c(ind, CR))
  
  # extract data from files and plot medians + interquartile ranges
  for (i in 1:ind) {
    for (cr in 1:CR) {
        Y_means[i, cr] <- mean(rel_Y[indices[i] - Time1 + 1, cr, ])
        SSB_means[i, cr] <- mean(rel_SSB[indices[i] - Time1 + 1, cr, ])
    }
  }
  
  data <- rep(2, 8)
  barplot(data, col = NA, border = FALSE, axes = TRUE)
  
  barplot(height = SSB_means, names.arg = 1:8, beside = TRUE, col = 'gray', 
          border = TRUE, xlab = 'Control Rule', ylim = c(0, 2))
  
  barplot(SSB_means[, 1])
  points(Y_means[, 1])
  

