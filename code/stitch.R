stitch <- function(Species, x, y) {
  
  library(abind)
  
  # put together filepaths for first set of sims
  file1 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_N.Rda', sep = '')
  file2 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_biomass.Rda', sep = '')
  # file3 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_SSB.Rda', sep = '')
  file4 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_yield.Rda', sep = '')
  # file5 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_effort.Rda', sep = '')
  # file6 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_DR.Rda', sep = '')
  # file7 <- paste('~/Projects/MS-thesis/data/', Species, '/', x, '_abundance.Rda', sep = '')
  
  # load them as objects
  load(file1); N1          <- sims_N  
  load(file2); biomass1    <- sims_biomass
  # load(file3); SSB1        <- sims_SSB
  load(file4); yield1      <- sims_yield
  # load(file5); effort1     <- sims_effort
  # load(file6); DR1         <- sims_DR
  # load(file7); abundance1  <- sims_abundance
  
  # put together filepaths for second set of sims
  file8 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_N.Rda', sep = '')
  file9 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_biomass.Rda', sep = '')
  # file10 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
  #                y, '_SSB.Rda', sep = '')
  file11 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                 y, '_yield.Rda', sep = '')
  # file12 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
  #                y, '_effort.Rda', sep = '')
  # file13 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
  #                y, '_DR.Rda', sep = '')
  # file14 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
  #                y, '_abundance.Rda', sep = '')
  
  # load them as objects
  load(file8);  N2          <- sims_N  
  load(file9);  biomass2    <- sims_biomass
  # load(file10); SSB2        <- sims_SSB
  load(file11); yield2      <- sims_yield
  # load(file12); effort2     <- sims_effort
  # load(file13); DR2         <- sims_DR
  # load(file14); abundance2  <- sims_abundance
  
  # combine objects into new ones
  sims_N         <- abind(N1, N2, along = 6)
  sims_biomass   <- abind(biomass1, biomass2, along = 5)
  # sims_SSB       <- abind(SSB1, SSB2, along = 5)
  sims_yield     <- abind(yield1, yield2, along = 5)
  # sims_effort    <- abind(effort1, effort2, along = 4)
  # sims_DR        <- abind(DR1, DR2, along = 4)
  # sims_abundance <- abind(abundance1, abundance2, along = 5)
  
  # new number of sims
  z = x + y
  
  # create new filepaths
  file15 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_N.Rda', sep = '')
  file16 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_biomass.Rda', sep = '')
  # file17 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
  #                 z, '_SSB.Rda', sep = '')
  file18 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
                  z, '_yield.Rda', sep = '')
  # file19 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
  #                 z, '_effort.Rda', sep = '')
  # file20 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
  #                 z, '_DR.Rda', sep = '')
  # file21 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/', Species, '/', 
  #                 z, '_abundance.Rda', sep = '')
  
  # save new objects to new filepaths
  save(sims_N, file = file15)
  save(sims_biomass, file = file16)
  # save(sims_SSB, file = file17)
  save(sims_yield, file = file18)
  # save(sims_effort, file = file19)
  # save(sims_DR, file = file20)
  # save(sims_abundance, file = file21)
  
  file.remove(file1)
  file.remove(file2)
  # file.remove(file3)
  file.remove(file4)
  # file.remove(file5)
  # file.remove(file6)
  # file.remove(file7)
  file.remove(file8)
  file.remove(file9)
  # file.remove(file10)
  file.remove(file11)
  # file.remove(file12)
  # file.remove(file13)
  # file.remove(file14)
  
}