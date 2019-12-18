selectivity_at_age <- function(fleets, L, max_age, rec_age, alpha, A50_up, 
                               A50_down, F_fin, beta, n, cf, age) {
  
  # length of fleet vector
  f <- length(fleets)
  
  # translate L50_up and L50_down to lengths instead of ages
  L50up <- L[A50_up - rec_age + 1]
  L50down <- L[A50_down - rec_age + 1]
  
  # initialize upcurves and downcurves
  upcurve <- array(rep(NA, f*(max_age + 1)), c(f, max_age + 1))
  downcurve <- array(rep(NA, f*(max_age + 1)), c(f, max_age + 1))
  
  for (i in 1:rec_age) {
    upcurve[, i] <- downcurve[, i] <- 0
  }
  
  # initialize selectivity at age array
  # dimensions = age * fleet
  S <- array(rep(0, f*n), c(f, n))
  
  for (i in 1:f) {
    upcurve[i, age + 1] <- 1 / (1 + exp(-1*alpha[i]*(L - L50up[i])))
    
    downcurve[i, age + 1] <- 1 - 
      (1 - F_fin[i]) / (1 + exp(-1*beta[i]*(L - L50down[i])))
    
    if (beta[i] == 0) {downcurve[i, age + 1] <- rep(1, n)}
    
    for (a in 1:n) {
      S[i, a] <- min(upcurve[i, a + rec_age], downcurve[i, a + rec_age])
    }
    
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
  # lines(age, S[2, ], type = 'l', lwd = 2, col = 'blue')
  # lines(age, S[3, ], type = 'l', lwd = 2, col = 'green')
  # # lines(age, S[4, ], type = 'l', lwd = 2, col = 'yellow')
  # # lines(age, S[5, ], type = 'l', lwd = 2, col = 'red')
  # legend(x = 'topright', fleets, lwd = 2, cex = 0.8,
  #        col = c('purple', 'blue', 'green')) #, 'yellow', 'red'))
  
  S[i, ] <- cf[i]*S[i, ]
  selectivity <- colSums(S)
  
  return (selectivity)
  
}