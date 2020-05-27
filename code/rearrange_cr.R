library(abind)
species_list = c('BR_OR_2015', 'CAB_OR_2019', 'LING_OW_2017', 'CR_OR_2015')

for (s in 1:length(species_list)) {
  
  Species = species_list[s]

  # put together filepaths for first set of sims
  file1 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/435_N.Rda', sep = '')
  file2 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/435_biomass.Rda', sep = '')
  file3 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/435_SSB.Rda', sep = '')
  file4 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/435_yield.Rda', sep = '')
  file5 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/435_effort.Rda', sep = '')
  file6 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/435_DR.Rda', sep = '')
  file7 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/435_abundance.Rda', sep = '')
  
  # load them as objects
  load(file1); N1          <- sims_N  
  load(file2); biomass1    <- sims_biomass
  load(file3); SSB1        <- sims_SSB
  load(file4); yield1      <- sims_yield
  load(file5); effort1     <- sims_effort
  load(file6); DR1         <- sims_DR
  load(file7); abundance1  <- sims_abundance
  
  # rearrange control rules
  sims_N         <- N1[, , , c(2, 5, 1, 4, 3, 6), , ]
  sims_biomass   <- biomass1[, , c(2, 5, 1, 4, 3, 6), , ]
  sims_SSB       <- SSB1[, , c(2, 5, 1, 4, 3, 6), , ]
  sims_yield     <- yield1[, c(2, 5, 1, 4, 3, 6), , ]
  sims_effort    <- effort1[, c(2, 5, 1, 4, 3, 6), , ]
  sims_DR        <- DR1[, c(2, 5, 1, 4, 3, 6), , ]
  sims_abundance <- abundance1[, , c(2, 5, 1, 4, 3, 6), , ]
  
  # create new filepaths
  file15 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_N.Rda', sep = '')
  file16 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_biomass.Rda', sep = '')
  file17 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_SSB.Rda', sep = '')
  file18 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_yield.Rda', sep = '')
  file19 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_effort.Rda', sep = '')
  file20 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_DR.Rda', sep = '')
  file21 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_abundance.Rda', sep = '')
  
  # save new objects to new filepaths
  save(sims_N, file = file15)
  save(sims_biomass, file = file16)
  save(sims_SSB, file = file17)
  save(sims_yield, file = file18)
  save(sims_effort, file = file19)
  save(sims_DR, file = file20)
  save(sims_abundance, file = file21)
  
  # file.remove(file1)
  # file.remove(file2)
  # file.remove(file3)
  # file.remove(file4)
  # file.remove(file5)
  # file.remove(file6)
  # file.remove(file7)
  
}
