# strip extra FDRs from Black Rockfish, Lingcod, Canary Rockfish

species_list = c('BR_OR_2015', 'LING_OW_2017', 'CR_OR_2015')

for (s in 1:length(species_list)) {
  
  Species <- species_list[s]
  
  load(paste("~/Projects/MS-thesis/data/Both/", Species, "/1166_abundance.Rda", sep = ''))
  load(paste("~/Projects/MS-thesis/data/Both/", Species, "/1166_biomass.Rda", sep = ''))
  load(paste("~/Projects/MS-thesis/data/Both/", Species, "/1166_DR.Rda", sep = ''))
  load(paste("~/Projects/MS-thesis/data/Both/", Species, "/1166_effort.Rda", sep = ''))
  load(paste("~/Projects/MS-thesis/data/Both/", Species, "/1166_N.Rda", sep = ''))
  load(paste("~/Projects/MS-thesis/data/Both/", Species, "/1166_SSB.Rda", sep = ''))
  load(paste("~/Projects/MS-thesis/data/Both/", Species, "/1166_yield.Rda", sep = ''))
  
  # load them as objects
  N1          <- sims_N  
  biomass1    <- sims_biomass
  SSB1        <- sims_SSB
  yield1      <- sims_yield
  effort1     <- sims_effort
  DR1         <- sims_DR
  abundance1  <- sims_abundance
  
  sims_N <- N1[, , , , 3:6, ]
  sims_biomass <- biomass1[, , , 3:6, ]
  sims_SSB <- SSB1[, , , 3:6, ]
  sims_yield <- yield1[, , 3:6, ]
  sims_effort <- effort1[, , 3:6, ]
  sims_DR <- DR1[, , 3:6, ]
  sims_abundance <- abundance1[, , , 3:6, ]
  
  # create new filepaths
  file1 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/1166_N.Rda', sep = '')
  file2 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/1166_biomass.Rda', sep = '')
  file3 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/1166_SSB.Rda', sep = '')
  file4 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/1166_yield.Rda', sep = '')
  file5 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/1166_effort.Rda', sep = '')
  file6 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/1166_DR.Rda', sep = '')
  file7 <- paste('~/Projects/MS-thesis/data/Both_before/', Species, '/1166_abundance.Rda', sep = '')
  
  # save new objects to new filepaths
  save(sims_N, file = file1)
  save(sims_biomass, file = file2)
  save(sims_SSB, file = file3)
  save(sims_yield, file = file4)
  save(sims_effort, file = file5)
  save(sims_DR, file = file6)
  save(sims_abundance, file = file7)
  
}

