replicate_base_model <- function(Species, num_sims) {
  
  library(remotes)
  remotes::install_github('vquennessen/densityratio')
  library(densityratio)
  library(abind)
  
  sims_N <- array()
  sims_biomass <- array()
  sims_SSB <- array()
  sims_yield <- array()
  sims_effort <- array()
  sims_DR <- array()
  
  run_base_model_once <- function(Species) {
    
    # Species = 'BR_OR_2015'
    # num_sims = 2
    
    output <- base_model(Species, R0 = 1e5, A = 5, MPA = 3, Time1 = 50, 
                         Time2 = 20, Recruitment_mode = 'pool', Error = 0.05, 
                         Stochasticity = FALSE, Surveys = TRUE, 
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
    
    
    sims_N <- ifelse(is.na(dim(sims_N)[2]), output$N, 
                     abind(sims_N, output$N, along = 6))
    sims_biomass <- ifelse(is.na(dim(sims_biomass)[2]), output$Biomass, 
                           abind(sims_biomass, output$Biomass, along = 5))
    sims_SSB <- ifelse(is.na(dim(sims_SSB)[2]), output$SSB, 
                       abind(sims_SSB, output$SSB, along = 5))
    sims_yield <- ifelse(is.na(dim(sims_yield)[2]), output$Yield, 
                         abind(sims_yield, output$Yield, along = 5))
    sims_effort <- ifelse(is.na(dim(sims_effort)[2]), output$Effort, 
                          abind(sims_effort, output$Effort, along = 4))
    sims_DR <- ifelse(is.na(dim(sims_DR)[2]), output$Density_ratio, 
                      abind(sims_DR, output$Density_ratio, along = 4))
    
  }
  
  replicate(num_sims, run_base_model_once(Species))
  
  Q <- ifelse(num_sims < 1000, num_sims,  paste("1e", log10(num_sims), sep = ''))
  
  filepath1 = paste('../data/', Species, '/', Q, '_N.Rda', sep = '')
  filepath2 = paste('../data/', Species, '/', Q, "_biomass.Rda", sep = '')
  filepath3 = paste('../data/', Species, '/', Q, "_SSB.Rda", sep = '')
  filepath4 = paste('../data/', Species, '/', Q, "_yield.Rda", sep = '')
  filepath5 = paste('../data/', Species, '/', Q, "_effort.Rda", sep = '')
  filepath6 = paste('../data/', Species, '/', Q, "_DR.Rda", sep = '')
  
  save(sims_N, file = filepath1)
  save(sims_biomass, file = filepath2)
  save(sims_SSB, file = filepath3)
  save(sims_yield, file = filepath4)
  save(sims_effort, file = filepath5)
  save(sims_DR, file = filepath6)
  
}


# timings:
source('run_base_model.R')
system.time(run_base_model(Species = 'BR_OR_2015', num_sims = 2))
system.time(lapply('BR_OR_2015', run_base_model, num_sims = 2))

source('replicate_base_model.R')
system.time(replicate_base_model(Species = 'BR_OR_2015', num_sims = 2))
