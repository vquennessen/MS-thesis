library(remotes)
remotes::install_github('vquennessen/densityratio')
library(densityratio)
source('./run_base_model.R')

Species = c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')
num_sims = 1e1

function2 <- function(Species, num_sims) {
  
  
  
  arrays <- lapply(Species, run_base_model)
  
  sims_N       <- arrays[[1]][[1]]
  sims_biomass <- arrays[[1]][[2]]
  sims_SSB     <- arrays[[1]][[3]]
  sims_yield   <- arrays[[1]][[4]]
  sims_effort  <- arrays[[1]][[5]]
  sims_DR      <- arrays[[1]][[6]]
  
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

lapply('BR_OR_2015', function2, num_sims)
