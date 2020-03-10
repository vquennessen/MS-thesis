cluster <- function(Species, Final_DR, num_sims) {
  
  # set arguments
  R0 = 1e5
  A = 5
  MPA = 3
  Time1 = 50
  Time2 = 20
  Error = 0.05
  Transects = 24
  Plotting = FALSE
  Control_rules = c(1:6)
  CR <- length(Control_rules)
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
  
  y_DR <- densityratio::transient_DR(Time1 = 50, TimeT = 70, Final_DR,
                                     Nat_mortality = c(0.09, 0.14, 0.19), nm = 2)
  Rec_age <- parameters(Species)[[3]]
  Max_age <- parameters(Species)[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  
  # total time
  TimeT <- Time1 + Time2
  
  # initialize yield and biomass arrays
  sims_N <- array(rep(0, n*MPA*TimeT*CR*1*num_sims), 
                  c(n, MPA, TimeT, CR, 1, num_sims))
  sims_biomass <- array(rep(0, MPA*TimeT*CR*1*num_sims), 
                        c(MPA, TimeT, CR, 1, num_sims))
  sims_SSB <- array(rep(0, MPA*TimeT*CR*1*num_sims), 
                    c(MPA, TimeT, CR, 1, num_sims))
  sims_yield <- array(rep(0, (MPA - 1)*TimeT*CR*num_sims), 
                      c((MPA - 1), TimeT, CR, num_sims))
  sims_Effort <- array(rep(0, (MPA - 1)*TimeT*CR*num_sims), 
                       c(MPA - 1, TimeT, CR, num_sims))
  sims_DR <- array(rep(0, (Time2 + 1)*CR*num_sims), c(Time2 + 1, CR, num_sims))
  
  # run the model for each simulation
  for (i in 1:num_sims) {
    
    output <- base_model(Species, R0, A, MPA, Time1, Time2, Recruitment_mode, 
                         Error, Stochasticity, Surveys, Fishery_management, 
                         Fishing, Transects, Adult_movement, Plotting, Final_DR, 
                         Years_sampled, Areas_sampled, Ind_sampled, Allocation, 
                         BM, Control_rules, Output.N = TRUE, Output.Biomass = TRUE, 
                         Output.SSB = TRUE, Output.Yield = TRUE, 
                         Output.Effort = TRUE, Output.Density.Ratio = TRUE)
    
    # save the relative yield and biomasses for all areas, times after reserve
    # implementation, and control rules
    sims_N[, , , , , i]     <- output$N
    sims_biomass[, , , , i] <- output$Biomass
    sims_SSB[, , , , i]     <- output$SSB
    sims_yield[, , , i]     <- output$Yield
    sims_Effort[, , , i]    <- output$Effort
    sims_DR[, , i]          <- output$Density_ratio
    
    if (i %% (num_sims/10) == 0) {
      percent <- 100*i/num_sims
      thing <- paste('Final DR ', Final_DR, ': ', percent, '% done.', sep = '')
      print(thing) 
    }  
  }
  
  Q <- ifelse(num_sims < 1000, num_sims,  paste("1e", log10(num_sims), sep = ''))
  
  filepath1 = paste('../data/', Species, '/E_', Q, "_", Final_DR, "_yield.Rda",
                    sep = '')
  filepath2 = paste('../data/', Species, '/E_', Q, "_", Final_DR,
                    "_biomass.Rda", sep = '')
  filepath3 = paste('../data/', Species, '/E_', Q, "_", Final_DR, "_SSB.Rda",
                    sep = '')
  filepath4 = paste('../data/', Species, '/E_', Q, "_", Final_DR, "_DR.Rda",
                    sep = '')
  filepath5 = paste('../data/', Species, '/E_', Q, '-', Final_DR, '_N.Rda',
                    sep = '')
  
  # filepath1 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/',
  #                   Species, '/ E_', Q, "_", Final_DR, "_yield.Rda", sep = '')
  # filepath2 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/',
  #                   Species, '/ E_', Q, "_", Final_DR, "_biomass.Rda", sep = '')
  # filepath3 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/',
  #                   Species, '/ E_', Q, "_", Final_DR, "_SSB.Rda", sep = '')
  # filepath4 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/',
  #                   Species, '/ E_', Q, "_", Final_DR, "_DR.Rda", sep = '')
  # filepath5 = paste('home/quennessenv/ExpanDrive/Box/Quennessen_Thesis/data/',
  #                   Species, '/ E_', Q, '-', Final_DR, '_N.Rda', sep = '')
  
  save(sims_yield, file = filepath1)
  save(sims_biomass, file = filepath2)
  save(sims_SSB, file = filepath3)
  save(sims_DR, file = filepath4)
  save(sims_N, file = filepath5)

}
