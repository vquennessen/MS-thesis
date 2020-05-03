var_stitch <- function(x, y) {
  
  library(abind)
  
  # put together filepaths for first set of sims
  file1 <- paste('~/Projects/MS-thesis/data/variance_runs/', x, '_biomass.Rda', sep = '')
  file2 <- paste('~/Projects/MS-thesis/data/variance_runs/', x, '_yield.Rda', sep = '')
  file3 <- paste('~/Projects/MS-thesis/data/variance_runs/', x, '_abundance.Rda', sep = '')
  
  # load them as objects
  load(file1); biomass1    <- sims_biomass
  load(file2); yield1      <- sims_yield
  load(file3); abundance1  <- sims_abundance
  
  # put together filepaths for second set of sims
  file4 <- paste('~/Projects/MS-thesis/data/variance_runs/', y, '_biomass.Rda', sep = '')
  file5 <- paste('~/Projects/MS-thesis/data/variance_runs/', y, '_yield.Rda', sep = '')
  file6 <- paste('~/Projects/MS-thesis/data/variance_runs/', y, '_abundance.Rda', sep = '')
  
  # load them as objects
  load(file4); biomass2    <- sims_biomass
  load(file5); yield2      <- sims_yield
  load(file6); abundance2  <- sims_abundance
  
  # combine objects into new ones
  sims_biomass   <- abind(biomass1, biomass2, along = 5)
  sims_yield     <- abind(yield1, yield2, along = 5)
  sims_abundance <- abind(abundance1, abundance2, along = 5)
  
  # new number of sims
  z = x + y
  
  # put together filepaths for second set of sims
  file7 <- paste('~/Projects/MS-thesis/data/variance_runs/', z, '_biomass.Rda', sep = '')
  file8 <- paste('~/Projects/MS-thesis/data/variance_runs/', z, '_yield.Rda', sep = '')
  file9 <- paste('~/Projects/MS-thesis/data/variance_runs/', z, '_abundance.Rda', sep = '')
  
  # save new objects to new filepaths
  save(sims_biomass, file = file7)
  save(sims_yield, file = file8)
  save(sims_abundance, file = file9)
  
  file.remove(file1)
  file.remove(file2)
  file.remove(file3)
  file.remove(file4)
  file.remove(file5)
  file.remove(file6)
  
}