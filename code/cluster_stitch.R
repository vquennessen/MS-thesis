cluster_stitch <- function(Scenario, x, y) {
  
  library(abind)
  
  if (Scenario == 'Variance') {
    species_list = c('CAB_OR_2019')
  } else {
    # species_list = c('CR_OR_2015', 'BR_OR_2015', 'CAB_OR_2019')
    species_list = c('LING_OW_2017', 'CR_OR_2015', 'BR_OR_2015', 'CAB_OR_2019')
  }
  
  for (s in 1:length(species_list)) {
    
    Species = species_list[s]
    
    # put together filepaths for first set of sims
    file1 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', x, '_biomass.Rda', sep = '')
    file2 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', x, '_yield.Rda', sep = '')
    
    if (Scenario != 'Variance') {
      file3 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', x, '_N.Rda', sep = '')
      file4 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', x, '_SSB.Rda', sep = '')
      file5 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', x, '_effort.Rda', sep = '')
      file6 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', x, '_DR.Rda', sep = '')
      file7 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', x, '_abundance.Rda', sep = '')
    }
    
    # load them as objects
    load(file1); biomass1    <- sims_biomass
    load(file2); yield1      <- sims_yield
    
    if (Scenario != 'Variance') {
      load(file3); N1          <- sims_N  
      load(file4); SSB1        <- sims_SSB
      load(file5); effort1     <- sims_effort
      load(file6); DR1         <- sims_DR
      load(file7); abundance1  <- sims_abundance      
    }
    
    # put together filepaths for second set of sims
    file8 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', y, '_biomass.Rda', sep = '')
    file9 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', y, '_yield.Rda', sep = '')
    
    if (Scenario != 'Variance') {
      file10 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', y, '_N.Rda', sep = '')
      file11 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', y, '_SSB.Rda', sep = '')
      file12 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', y, '_effort.Rda', sep = '')
      file13 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', y, '_DR.Rda', sep = '')
      file14 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', y, '_abundance.Rda', sep = '')
    }
    
    # load them as objects
    load(file8);  biomass2    <- sims_biomass
    load(file9); yield2      <- sims_yield
    
    if (Scenario != 'Variance') {
      load(file10);  N2          <- sims_N  
      load(file11); SSB2        <- sims_SSB
      load(file12); effort2     <- sims_effort
      load(file13); DR2         <- sims_DR
      load(file14); abundance2  <- sims_abundance      
    }
    
    # combine objects into new ones
    sims_biomass   <- abind(biomass1, biomass2, along = 5)
    sims_yield     <- abind(yield1, yield2, along = 4)
    
    if (Scenario != 'Variance') {
      sims_N         <- abind(N1, N2, along = 6)
      sims_SSB       <- abind(SSB1, SSB2, along = 5)
      sims_effort    <- abind(effort1, effort2, along = 4)
      sims_DR        <- abind(DR1, DR2, along = 4)
      sims_abundance <- abind(abundance1, abundance2, along = 5)      
    }
    
    # new number of sims
    z = x + y
    
    # create new filepaths
    file15 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', z, '_biomass.Rda', sep = '')
    file16 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', z, '_yield.Rda', sep = '')
    
    if (Scenario != 'Variance') {
      file17 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', z, '_N.Rda', sep = '')
      file18 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', z, '_SSB.Rda', sep = '')
      file19 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', z, '_effort.Rda', sep = '')
      file20 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', z, '_DR.Rda', sep = '')
      file21 <- paste('~/Documents/MS-thesis/data/', Scenario, '/', Species, '/', z, '_abundance.Rda', sep = '')      
    }
    
    # save new objects to new filepaths
    save(sims_biomass, file = file15)
    save(sims_yield, file = file16)
    
    if (Scenario != 'Variance') {
      save(sims_N, file = file17)
      save(sims_SSB, file = file18)
      save(sims_effort, file = file19)
      save(sims_DR, file = file20)
      save(sims_abundance, file = file21)      
    }
    
    # remove old files
    file.remove(file1)
    file.remove(file2)
    file.remove(file8)
    file.remove(file9)    
    
    if (Scenario != 'Variance'){
      file.remove(file3)
      file.remove(file4)
      file.remove(file5)
      file.remove(file6)
      file.remove(file7)      
      file.remove(file10)
      file.remove(file11)
      file.remove(file12)
      file.remove(file13)
      file.remove(file14)    
    }
    
  }
  
}
