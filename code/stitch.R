stitch <- function(Species, x, y) {
  
  library(abind)
  
  # put together filepaths for first set of sims
  file1 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_N.Rda', sep = '')
  file2 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_biomass.Rda', sep = '')
  file3 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_SSB.Rda', sep = '')
  file4 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_yield.Rda', sep = '')
  file5 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_effort.Rda', sep = '')
  file6 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_DR.Rda', sep = '')
  
  # load them as objects
  load(file1); N1       <- sims_N  
  load(file2); biomass1 <- sims_biomass
  load(file3); SSB1     <- sims_SSB
  load(file4); yield1   <- sims_yield
  load(file5); effort1  <- sims_effort
  load(file6); DR1      <- sims_DR
  
  # put together filepaths for second set of sims
  file7 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_N.Rda', sep = '')
  file8 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_biomass.Rda', sep = '')
  file9 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_SSB.Rda', sep = '')
  file10 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_yield.Rda', sep = '')
  file11 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_effort.Rda', sep = '')
  file12 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_DR.Rda', sep = '')
  
  # load them as objects
  load(file7);  N2       <- sims_N  
  load(file8);  biomass2 <- sims_biomass
  load(file9);  SSB2     <- sims_SSB
  load(file10); yield2   <- sims_yield
  load(file11); effort2  <- sims_effort
  load(file12); DR2      <- sims_DR
  
  # combine objects into new ones
  N       <- abind(N1, N2, along = 6)
  biomass <- abind(biomass1, biomass2, along = 5)
  SSB     <- abind(SSB1, SSB2, along = 5)
  yield   <- abind(yield1, yield2, along = 5)
  effort  <- abind(effort1, effort2, along = 4)
  DR      <- abind(DR1, DR2, along = 4)
  
  # new number of sims
  # z = x + y
  
  # create new filepaths
  file13 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_N.Rda', sep = '')
  file14 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_biomass.Rda', sep = '')
  file15 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_SSB.Rda', sep = '')
  file16 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_yield.Rda', sep = '')
  file17 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_effort.Rda', sep = '')
  file18 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_DR.Rda', sep = '')
  
  # save new objects to new filepaths
  save(N, file = file13)
  save(biomass, file = file14)
  save(SSB, file = file15)
  save(yield, file = file16)
  save(effort, file = file17)
  save(DR, file = file18)
  
}