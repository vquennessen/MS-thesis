selectivity_at_age <- function(Fleets, L, Max_age, Rec_age, Alpha, A50_up, 
                               A50_down, F_fin, Beta, Num, Cf, Age) {
  
  # length of fleet vector
  f <- length(Fleets)
  
  # translate L50_up and L50_down to lengths instead of ages
  L50up <- L[A50_up - Rec_age + 1]
  L50down <- L[A50_down - Rec_age + 1]
  
  # initialize upcurves and downcurves
  upcurve <- array(rep(NA, f*(Max_age + 1)), c(f, Max_age + 1))
  downcurve <- array(rep(NA, f*(Max_age + 1)), c(f, Max_age + 1))
  
  for (i in 1:Rec_age) {
    upcurve[, i] <- downcurve[, i] <- 0
  }
  
  # initialize selectivity at age array
  # dimensions = age * fleet
  S <- array(rep(0, f*Num), c(f, Num))
  
  for (i in 1:f) {
    upcurve[i, Age + 1] <- 1 / (1 + exp(-1*Alpha[i]*(L - L50up[i])))
    
    downcurve[i, Age + 1] <- 1 - 
      (1 - F_fin[i]) / (1 + exp(-1*Beta[i]*(L - L50down[i])))
    
    if (Beta[i] == 0) { downcurve[i, Age + 1] <- rep(1, Num) }
    
    for (a in 1:Num) {
      S[i, a] <- min(upcurve[i, a + Rec_age], downcurve[i, a + Rec_age])
    }
    
    S[i, ] <- Cf[i]*S[i, ]
    
  }
  
  # ############################################################################
  # ##### plot selectivities to double check they're right
  # ############################################################################
  # 
  # plot(age, S[1, ], type = 'l', lwd = 2, col = 'purple',
  #      ylim = c(0, 1),
  #      main = "Vic's Attempt",
  #      xlab = 'Age (year)',
  #      ylab = 'Selectivity',
  #      xlim = c(0, 40))
  # lines(Age, S[2, ], type = 'l', lwd = 2, col = 'blue')
  # lines(Age, S[3, ], type = 'l', lwd = 2, col = 'green')
  # # lines(Age, S[4, ], type = 'l', lwd = 2, col = 'yellow')
  # # lines(Age, S[5, ], type = 'l', lwd = 2, col = 'red')
  # legend(x = 'topright', Fleets, lwd = 2, cex = 0.8,
  #        col = c('purple', 'blue', 'green')) #, 'yellow', 'red'))
  
  selectivity <- colSums(S)
  
  return (selectivity)
  
}