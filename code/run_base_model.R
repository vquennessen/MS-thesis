rm(list = ls())
setwd('~/Projects/MS-thesis/code')
source('./plot_stuff.R')
devtools::install_git('https://github.com/vquennessen/densityratio.git')
library(densityratio)
  
  # set numbers of simulations
  num_sims <- 1
  
  # set arguments
  Species = 'BR_OR_2015'
  R0 = 1e5
  A = 5
  MPA = 3
  Time1 = 50
  Time2 = 20
  Error = 0.05
  Transects = 24
  Plotting = FALSE
  Final_DRs = c(0.2, 0.6, 1.0)
  Control_rules = c(1:6)
  Plot_individual_runs <- F
  BM <- FALSE
  Output.Yield = TRUE
  Output.Biomass = TRUE
  Output.SSB = TRUE
  Output.N = TRUE
  Output.Density.Ratio = TRUE
  Output.Effort = TRUE
  Surveys = TRUE
  Fishing = TRUE
  Stochasticity = TRUE
  Recruitment_mode = 'pool'
  Allocation = 'IFD'
  Adult_movement = TRUE
  Fishery_management = TRUE
  Years_sampled = 1
  Areas_sampled = 'all'
  Ind_sampled = 'all'
  
  # dimensions
  TimeT <- Time1 + Time2
  CR <- length(Control_rules)
  FDR <- length(Final_DRs)
  Rec_age <- parameters(Species)[[3]]
  Max_age <- parameters(Species)[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  
  # initialize yield and biomass arrays
  sims_N <- array(rep(0, n*MPA*TimeT*CR*FDR*num_sims), 
                  c(n, MPA, TimeT, CR, FDR, num_sims))
  sims_biomass <- array(rep(0, MPA*TimeT*CR*FDR*num_sims), 
                        c(MPA, TimeT, CR, FDR, num_sims))
  sims_SSB <- array(rep(0, MPA*TimeT*CR*FDR*num_sims), 
                    c(MPA, TimeT, CR, FDR, num_sims))
  sims_yield <- array(rep(0, (MPA - 1)*TimeT*CR*FDR*num_sims), 
                      c((MPA - 1), TimeT, CR, FDR, num_sims))
  sims_effort <- array(rep(0, (MPA - 1)*TimeT*CR*FDR*num_sims), 
                       c(MPA - 1, TimeT, CR, FDR, num_sims))
  sims_DR <- array(rep(0, (Time2 + 1)*CR*FDR*num_sims), 
                   c(Time2 + 1, CR, FDR, num_sims))
  
  # run the model for each simulation
  system.time(
    for (i in 1:num_sims) {
      
      output <- base_model(Species, R0, A, MPA, Time1, Time2, Recruitment_mode, 
                           Error, Stochasticity, Surveys, Fishery_management, 
                           Fishing, Transects, Adult_movement, Plotting, Final_DRs, 
                           Years_sampled, Areas_sampled, Ind_sampled, Allocation, 
                           BM, Control_rules, Output.N = TRUE, Output.Biomass = TRUE, 
                           Output.SSB = TRUE, Output.Yield = TRUE, 
                           Output.Effort = TRUE, Output.Density.Ratio = TRUE)
      
      # save the relative yield and biomasses for all areas, times after reserve
      # implementation, and control rules
      sims_N[, , , , , i]       <- output$N
      sims_biomass[, , , , i]   <- output$Biomass
      sims_SSB[, , , , i]       <- output$SSB
      sims_yield[, , , , i]     <- output$Yield
      sims_effort[, , , , i]    <- output$Effort
      sims_DR[, , , i]          <- output$Density_ratio
      
      if (i %% (num_sims/1) == 0) {
        percent <- 100*i/num_sims
        thing <- paste(Species, ': ', percent, '% done.', sep = '')
        print(thing) 
      }  
    }
  )
  
  Q <- ifelse(num_sims < 1000, num_sims,  paste("1e", log10(num_sims), sep = ''))
  
  filepath1 = paste('../data/', Species, '/E_', Q, "_yield.Rda", sep = '')
  filepath2 = paste('../data/', Species, '/E_', Q, "_biomass.Rda", sep = '')
  filepath3 = paste('../data/', Species, '/E_', Q, "_SSB.Rda", sep = '')
  filepath4 = paste('../data/', Species, '/E_', Q, "_DR.Rda", sep = '')
  filepath5 = paste('../data/', Species, '/E_', Q, '_N.Rda', sep = '')
  
  # filepath1 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/', 
  #                   Species, '/E_', Q, "_yield.Rda", sep = '')
  # filepath2 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/', 
  #                   Species, '/E_', Q, "_biomass.Rda", sep = '')
  # filepath3 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/', 
  #                   Species, '/E_', Q, "_SSB.Rda", sep = '')
  # filepath4 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/', 
  #                   Species, '/E_', Q, "_DR.Rda", sep = '')
  # filepath5 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/', 
  #                   Species, '/E_', Q, '_N.Rda', sep = '')
  
  save(sims_yield, file = filepath1)
  save(sims_biomass, file = filepath2)
  save(sims_SSB, file = filepath3)
  save(sims_DR, file = filepath4)
  save(sims_N, file = filepath5)
