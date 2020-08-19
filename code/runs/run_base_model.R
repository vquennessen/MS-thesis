run_base_model <- function(Species, num_sims, Scenario, Final_DRs) {
  
  library(remotes)
  remotes::install_github('vquennessen/densityratio')
  library(densityratio)
  
  start_time <- paste('Start time: ', Sys.time(), ' - ', Species, ' - ',
                      num_sims, ' sims', sep = '')
  write(start_time, file = 'progress.txt', append = TRUE)
  
  if (Scenario == 'None' | Scenario == 'Test' | Scenario == 'Areas') {
    Sampling_Var = FALSE
    Transects = 24
    Recruitment_Var = FALSE
    
  } else if (Scenario == 'Sampling.LowT') {
    Sampling_Var = TRUE
    Transects = 12
    Recruitment_Var = FALSE
    
  } else if (Scenario == 'New.Sampling') {
    Sampling_Var = TRUE
    Transects = 24
    Recruitment_Var = FALSE
    
  } else if (Scenario == 'Sampling.HighT') {
    Sampling_Var = TRUE
    Transects = 36
    Recruitment_Var = FALSE
    
  } else if (Scenario == 'Recruitment') {
    Sampling_Var = FALSE
    Transects = 24
    Recruitment_Var = TRUE
    
  } else if (Scenario == 'Both.LowT') {
    Sampling_Var = TRUE
    Transects = 12
    Recruitment_Var = TRUE
    
  } else if (Scenario == 'Both' | Scenario == 'Variance') {
    Sampling_Var = TRUE
    Transects = 24
    Recruitment_Var = TRUE
    
  } else if (Scenario == 'Both.HighT') {
    Sampling_Var = TRUE
    Transects = 36
    Recruitment_Var = TRUE
  } 
  
  # set arguments
  if (Scenario != 'Variance') {
    Output.N = TRUE
    Output.Abundance = TRUE
    Output.Biomass = TRUE
    Output.SSB = TRUE
    Output.Yield = TRUE
    Output.Effort = TRUE
    Output.Density.Ratio = TRUE
  } else {
    Output.N = FALSE
    Output.Abundance = FALSE
    Output.Biomass = TRUE
    Output.SSB = FALSE
    Output.Yield = TRUE
    Output.Effort = FALSE
    Output.Density.Ratio = FALSE
  }
  
  R0 = 1e+5
  A = 5
  MPA = 3
  Time1 = 50
  Time2 = 20
  TimeT = Time1 + Time2
  Recruitment_mode = 'pool'
  Surveys = TRUE
  Fishery_management = TRUE
  Fishing = TRUE
  Adult_movement = TRUE
  Years_sampled = 1
  Areas_sampled = 'all'
  Ind_sampled = 'all'
  Floor_DR = 0.2
  Allocation = 'IFD'
  BM = FALSE
  LDP = 0.1
  M_Error = 0.05
  
  # dimensions
  CR <- 6
  NM <- 2
  FDR <- length(Final_DRs)
  Rec_age <- parameters(Species)[[3]]
  Max_age <- parameters(Species)[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  
  # initialize yield and biomass arrays
  sims_biomass <- array(rep(0, A*(Time2 + 1)*CR*FDR*num_sims),
                        c(A, Time2 + 1, CR, FDR, num_sims))
  sims_yield <- array(rep(0, (Time2 + 1)*CR*FDR*num_sims),
                      c(Time2 + 1, CR, FDR, num_sims))
  
  if (Scenario != 'Variance') {
    sims_N <- array(rep(0, n*A*(Time2 + 1)*CR*FDR*num_sims),
                    c(n, A, Time2 + 1, CR, FDR, num_sims))
    sims_DR <- array(rep(0, (Time2 + 1)*CR*FDR*num_sims),
                     c(Time2 + 1, CR, FDR, num_sims))
    sims_SSB <- array(rep(0, A*(Time2 + 1)*CR*FDR*num_sims),
                      c(A, Time2 + 1, CR, FDR, num_sims))
    sims_effort <- array(rep(0, TimeT*CR*FDR*num_sims),
                         c(TimeT, CR, FDR, num_sims))
    sims_abundance <- array(rep(0, A*(Time2 + 1)*CR*FDR*num_sims),
                            c(A, Time2 + 1, CR, FDR, num_sims))
  }
  
  # run the model for each simulation
  for (i in 1:num_sims) {
    
    output <- base_model(Species, R0, A, MPA, Time1, Time2, Recruitment_mode, 
                         M_Error, Sampling_Var, Recruitment_Var, Surveys, 
                         Fishery_management, Fishing, Transects, Adult_movement, 
                         Final_DRs, Years_sampled, Areas_sampled, Ind_sampled, 
                         Floor_DR, Allocation, BM, LDP, Output.N, 
                         Output.Abundance, Output.Biomass, Output.SSB, 
                         Output.Yield, Output.Effort, Output.Density.Ratio)
    
    # save the relative yield and biomasses for all areas, times after reserve
    # implementation, control rules, and final density ratios
    sims_biomass[, , , , i]   <- output$Biomass
    sims_yield[, , , i]       <- output$Yield
    
    if (Scenario != 'Variance') {
      sims_N[, , , , , i]       <- output$N
      sims_DR[, , , i]          <- output$Density_ratio
      sims_SSB[, , , , i]       <- output$SSB
      sims_effort[, , , i]      <- output$Effort
      sims_abundance[, , , , i] <- output$Abundance
    }
    
    if (i %% (num_sims/10) == 0) {
      update <- paste(Sys.time(), ' - ', Species, ' - ', i/num_sims*100, 
                      '% done!', sep = '')
      write(update, file = 'progress.txt', append = TRUE)
    }
    
  }
  
  # get filepaths to save objects to
  # Q <- ifelse(num_sims < 1000, num_sims,  paste("1e", log10(num_sims), sep = ''))
  Q <- num_sims
  
  filepath1 = paste('../data/', Scenario, '/', Species, '/', Q, '_biomass.Rda', sep = '')
  filepath2 = paste('../data/', Scenario, '/', Species, '/', Q, '_yield.Rda', sep = '')
  
  if (Scenario != 'Variance') {
    filepath3 = paste('../data/', Scenario, '/', Species, '/', Q, '_N.Rda', sep = '')
    filepath4 = paste('../data/', Scenario, '/', Species, '/', Q, '_DR.Rda', sep = '')
    filepath5 = paste('../data/', Scenario, '/', Species, '/', Q, '_effort.Rda', sep = '')
    filepath6 = paste('../data/', Scenario, '/', Species, '/', Q, '_SSB.Rda', sep = '')
    filepath7 = paste('../data/', Scenario, '/', Species, '/', Q, '_abundance.Rda', sep = '')
  }
  
  # save objects
  save(sims_biomass, file = filepath1)
  save(sims_yield, file = filepath2)
  
  if (Scenario != 'Variance') {
    save(sims_N, file = filepath3)
    save(sims_DR, file = filepath4)
    save(sims_effort, file = filepath5)
    save(sims_SSB, file = filepath6)
    save(sims_abundance, file = filepath7)
  }
  
}
