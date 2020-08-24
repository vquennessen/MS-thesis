version2 <- function(Scenario, x) {
  
  if (Scenario == 'Variance') {
    species_list = c('CAB_OR_2019')
  } else {
    # species_list = c('LING_OW_2017', 'CR_OR_2015', 'BR_OR_2015', 'CAB_OR_2019')
    species_list = c('BR_OR_2015')
    
  }
  
  for (s in 1:length(species_list)) {
    
    Species = species_list[s]
    
    # put together filepaths for first set of sims
    file1 <- paste('~/Projects/MS-thesis/data/', Scenario, '/', Species, '/', x, '_biomass.Rda', sep = '')
    file2 <- paste('~/Projects/MS-thesis/data/', Scenario, '/', Species, '/', x, '_yield.Rda', sep = '')
    
    if (Scenario != 'Variance') {
      file3 <- paste('~/Projects/MS-thesis/data/', Scenario, '/', Species, '/', x, '_N.Rda', sep = '')
      file4 <- paste('~/Projects/MS-thesis/data/', Scenario, '/', Species, '/', x, '_SSB.Rda', sep = '')
      file5 <- paste('~/Projects/MS-thesis/data/', Scenario, '/', Species, '/', x, '_effort.Rda', sep = '')
      file6 <- paste('~/Projects/MS-thesis/data/', Scenario, '/', Species, '/', x, '_DR.Rda', sep = '')
      file7 <- paste('~/Projects/MS-thesis/data/', Scenario, '/', Species, '/', x, '_abundance.Rda', sep = '')
    }
    
    # load them as objects
    load(file1)
    load(file2)
    
    if (Scenario != 'Variance') {
      load(file3)
      load(file4)
      load(file5)
      load(file6)
      load(file7)    
    }
    
    # save new objects to new filepaths
    save(sims_biomass, file = file1, version = 2)
    save(sims_yield, file = file2, version = 2)
    
    if (Scenario != 'Variance') {
      save(sims_N, file = file3, version = 2)
      save(sims_SSB, file = file4, version = 2)
      save(sims_effort, file = file5, version = 2)
      save(sims_DR, file = file6, version = 2)
      save(sims_abundance, file = file7, version = 2)      
    }
    
  }
  
}