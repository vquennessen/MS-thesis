plot_stuff <- function(filepath1, filepath3, A, Time1, Time2, CR, num_sims, 
                       sample_size, PD, plot_individual_runs, y_DR, species, 
                       final_DR) {
  
  # load files
  load(filepath1)
  load(filepath3)
  
  # starting yield and biomass
  start_Y <- array(rep(0, CR*num_sims), c(CR, num_sims))
  start_SSB <- array(rep(0, CR*num_sims), c(CR, num_sims))
  
  for (cr in 1:CR) {
    for (sim in 1:num_sims) {
        start_Y <- sum(sims_yield[, Time1, cr, sim])
        start_SSB <- sum(sims_SSB[, Time1, cr, sim])
    }
  }

  # pull out only every 5 years
  TimeT <- Time1 + Time2
  indices <- seq(Time1 + 1, TimeT, by = 5)
  ind <- length(indices)
  
  # save objects
  Y <- sims_yield[, indices, , ]
  SSB <- sims_SSB[, indices, , ]

  # replace NA's in sims_yield with 0s
  Y[is.na(Y)] <- 0
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- array(rep(NA, ind*CR), c(ind, CR))
  SSB_medians <- array(rep(NA, ind*CR), c(ind, CR))

  # extract data from files and plot medians + interquartile ranges
  for (t in 1:ind) {
    for (cr in 1:CR) {
      for (a in 1:A) {
        
        Y_medians[t, cr] <- median(Y[, t, cr, ])
        SSB_medians[a, t, cr] <- median(SSB[a, t, cr, ])
  
      }
    }
  }
  
  barplot(SSB_medians)
  
  
  
  }

  