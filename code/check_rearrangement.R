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
  
  # create new filepaths
  file15 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_N.Rda', sep = '')
  file16 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_biomass.Rda', sep = '')
  file17 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_SSB.Rda', sep = '')
  file18 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_yield.Rda', sep = '')
  file19 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_effort.Rda', sep = '')
  file20 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_DR.Rda', sep = '')
  file21 <- paste('C:/Users/Vic/Documents/Projects/MS-thesis/data/Both/', Species, '/435_abundance.Rda', sep = '')
  
  # load them as objects
  load(file15); N2          <- sims_N  
  load(file16); biomass2    <- sims_biomass
  load(file17); SSB2        <- sims_SSB
  load(file18); yield2      <- sims_yield
  load(file19); effort2     <- sims_effort
  load(file20); DR2         <- sims_DR
  load(file21); abundance2  <- sims_abundance
  
  # print differences
  print(Species)
  A <- c(sum(N2[, , , 1, , ] - N1[, , , 2, , ]), 
         sum(N2[, , , 2, , ] - N1[, , , 5, , ]), 
         sum(N2[, , , 3, , ] - N1[, , , 1, , ]),
         sum(N2[, , , 4, , ] - N1[, , , 4, , ]),
         sum(N2[, , , 5, , ] - N1[, , , 3, , ]),
         sum(N2[, , , 6, , ] - N1[, , , 6, , ]))
  A <- any(A != 0)
  if (A == 0) {print('N: all good') } else {print('N: Problem')}
  
  B <- c(sum(biomass2[, , 1, , ] - biomass1[, , 2, , ]), 
         sum(biomass2[, , 2, , ] - biomass1[, , 5, , ]), 
         sum(biomass2[, , 3, , ] - biomass1[, , 1, , ]),
         sum(biomass2[, , 4, , ] - biomass1[, , 4, , ]),
         sum(biomass2[, , 5, , ] - biomass1[, , 3, , ]),
         sum(biomass2[, , 6, , ] - biomass1[, , 6, , ]))
  B <- any(B != 0)
  if (B == 0) {print('Biomass: all good') } else {print('Biomass: Problem')}

  C <- c(sum(SSB2[, , 1, , ] - SSB1[, , 2, , ]), 
         sum(SSB2[, , 2, , ] - SSB1[, , 5, , ]), 
         sum(SSB2[, , 3, , ] - SSB1[, , 1, , ]),
         sum(SSB2[, , 4, , ] - SSB1[, , 4, , ]),
         sum(SSB2[, , 5, , ] - SSB1[, , 3, , ]),
         sum(SSB2[, , 6, , ] - SSB1[, , 6, , ]))
  C <- any(C != 0)
  if (C == 0) {print('SSB: all good') } else {print('SSB: Problem')}

  D <- c(sum(yield2[, 1, , ] - yield1[, 2, , ]), 
         sum(yield2[, 2, , ] - yield1[, 5, , ]), 
         sum(yield2[, 3, , ] - yield1[, 1, , ]),
         sum(yield2[, 4, , ] - yield1[, 4, , ]),
         sum(yield2[, 5, , ] - yield1[, 3, , ]),
         sum(yield2[, 6, , ] - yield1[, 6, , ]))
  D <- any(D != 0)
  if (D == 0) {print('Yield: all good') } else {print('Yield: Problem')}

  E <- c(sum(effort2[, 1, , ] - effort1[, 2, , ]), 
         sum(effort2[, 2, , ] - effort1[, 5, , ]), 
         sum(effort2[, 3, , ] - effort1[, 1, , ]),
         sum(effort2[, 4, , ] - effort1[, 4, , ]),
         sum(effort2[, 5, , ] - effort1[, 3, , ]),
         sum(effort2[, 6, , ] - effort1[, 6, , ]))
  E <- any(E != 0)
  if (E == 0) {print('Effort: all good') } else {print('Effort: Problem')}

  G <- c(sum(DR2[, 1, , ] - DR1[, 2, , ]), 
         sum(DR2[, 2, , ] - DR1[, 5, , ]), 
         sum(DR2[, 3, , ] - DR1[, 1, , ]),
         sum(DR2[, 4, , ] - DR1[, 4, , ]),
         sum(DR2[, 5, , ] - DR1[, 3, , ]),
         sum(DR2[, 6, , ] - DR1[, 6, , ]))
  G <- any(G != 0)
  if (G == 0) {print('DR: all good') } else {print('DR: Problem')}

  H <- c(sum(abundance2[, , 1, , ] - abundance1[, , 2, , ]), 
         sum(abundance2[, , 2, , ] - abundance1[, , 5, , ]), 
         sum(abundance2[, , 3, , ] - abundance1[, , 1, , ]),
         sum(abundance2[, , 4, , ] - abundance1[, , 4, , ]),
         sum(abundance2[, , 5, , ] - abundance1[, , 3, , ]),
         sum(abundance2[, , 6, , ] - abundance1[, , 6, , ]))
  H<- any(H != 0)
  if (H == 0) {print('Abundance: all good') } else {print('Abundance: Problem')}

}
