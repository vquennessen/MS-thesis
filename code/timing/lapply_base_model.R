lapply_base_model <- function(Species, num_sims) {
  
  library(remotes)
  remotes::install_github('vquennessen/densityratio')
  library(densityratio)
  
  A = 5
  MPA = 3
  Time1 = 50
  Time2 = 20
  Final_DRs = c(0.2, 0.4, 0.6, 0.8, 1)
  Control_rules = c(1:6)
  
  # dimensions
  TimeT <- Time1 + Time2
  CR <- length(Control_rules)
  FDR <- length(Final_DRs)
  Rec_age <- parameters(Species)[[3]]
  Max_age <- parameters(Species)[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  
  sims_N <- array(rep(0, n*MPA*TimeT*CR*FDR*num_sims),
                  c(n, MPA, TimeT, CR, FDR, num_sims))
  sims_biomass <- array(rep(0, MPA*TimeT*CR*FDR*num_sims),
                        c(MPA, TimeT, CR, FDR, num_sims))
  sims_SSB <- array(rep(0, MPA*TimeT*CR*FDR*num_sims),
                    c(MPA, TimeT, CR, FDR, num_sims))
  sims_yield <- array(rep(0, (MPA - 1)*TimeT*CR*FDR*num_sims),
                      c((MPA - 1), TimeT, CR, FDR, num_sims))
  sims_effort <- array(rep(0, TimeT*CR*FDR*num_sims),
                       c(TimeT, CR, FDR, num_sims))
  sims_DR <- array(rep(0, TimeT*CR*FDR*num_sims),
                   c(TimeT, CR, FDR, num_sims))
  
  run_base_model_once <- function(i, Species) { #}, sims_biomass, sims_SSB, 
                                  # sims_yield, sims_effort, sims_DR) {
    
    output <- base_model(Species, R0 = 1e5, A = 5, MPA = 3, Time1 = 50, 
                         Time2 = 20, Recruitment_mode = 'pool', Error = 0.05, 
                         Stochasticity = TRUE, Surveys = TRUE, 
                         Fishery_management = TRUE, Fishing = TRUE, 
                         Transects = 24, Adult_movement = TRUE, Plotting = FALSE, 
                         Final_DRs = c(0.2, 0.4, 0.6, 0.8, 1), Years_sampled = 1, 
                         Areas_sampled = 'all', Ind_sampled = 'all', 
                         Allocation = 'IFD', BM = FALSE, LDP = 0.1, 
                         Control_rules = c(1:6), Output.FM = FALSE, 
                         Output.N = TRUE, Output.Abundance = FALSE, 
                         Output.Biomass = TRUE, Output.SSB = TRUE, 
                         Output.Catch = FALSE, Output.Yield = TRUE, 
                         Output.Effort = TRUE, Output.Density.Ratio = TRUE)
    
    # save the relative yield and biomasses for all areas, times after reserve
    # implementation, control rules, and final density ratios
    sims_N[, , , , , i]       <- output$N
    sims_biomass[, , , , i]   <- output$Biomass
    sims_SSB[, , , , i]       <- output$SSB
    sims_yield[, , , , i]     <- output$Yield
    sims_effort[, , , i]      <- output$Effort
    sims_DR[, , , i]          <- output$Density_ratio
    
    if (i %% (num_sims/100) == 0) {
      update <- paste(Sys.time(), ' - ', Species, ' - ', i/num_sims*100, 
                      '% done!', sep = '')
      write(update, file = 'progress.txt', append = TRUE)
    }
    
    return(list(sims_N, sims_biomass, sims_SSB, sims_yield, sims_effort, 
                sims_DR))
    
  }
  
  returned <- lapply(1:num_sims, run_base_model_once, Species)
  
  sims_N       <- returned[[1]][[1]]
  sims_biomass <- returned[[1]][[2]]
  sims_SSB     <- returned[[1]][[3]]
  sims_yield   <- returned[[1]][[4]]
  sims_effort  <- returned[[1]][[5]]
  sims_DR      <- returned[[1]][[6]]
  
  # extract first number of num_sims and put num_sims label together
  z <- substr(num_sims, start = 0, stop = 1)
  Q <- ifelse(num_sims < 1000, num_sims,  paste(z, "e", log10(num_sims), sep = ''))
  
  # use num_sims label to set filepaths for different objects
  filepath1 = paste('../data/', Species, '/', Q, '_N.Rda', sep = '')
  filepath2 = paste('../data/', Species, '/', Q, "_biomass.Rda", sep = '')
  filepath3 = paste('../data/', Species, '/', Q, "_SSB.Rda", sep = '')
  filepath4 = paste('../data/', Species, '/', Q, "_yield.Rda", sep = '')
  filepath5 = paste('../data/', Species, '/', Q, "_effort.Rda", sep = '')
  filepath6 = paste('../data/', Species, '/', Q, "_DR.Rda", sep = '')
  
  # save objects to filepaths
  save(sims_N, file = filepath1)
  save(sims_biomass, file = filepath2)
  save(sims_SSB, file = filepath3)
  save(sims_yield, file = filepath4)
  save(sims_effort, file = filepath5)
  save(sims_DR, file = filepath6)
  
}
