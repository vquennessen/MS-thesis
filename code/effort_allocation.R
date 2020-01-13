effort_allocation <- function(t, cr, nm, Allocation, E, Yield, Time1, Inside, 
                              Outside) {
  
  #' number of areas not in a reserve
  out <- length(Outside)
  ins <- length(Inside)
  all <- out + ins
  
  #' If effort is allocated using the ideal free distribution, effort for one 
  #' year depends on the distribution of yield from the previous year
  if (Allocation == 'IFD') {
      
      prop_yield_outside <- Yield[Outside, t - 1, cr, nm]/sum(Yield[Outside, t - 1, cr, nm])
      
      E[Outside, t, cr, nm] <- sum(E[, t - 1, cr, nm])*prop_yield_outside
      E[Inside, t, cr, nm] <- 0
    
  #' Otherwise, distribute effort equally between the four areas outside the 
  #' marine reserve, regardless of yield
  } else if (Allocation == 'equal') {
    
    if (t < Time1) {
      
      E[, t, cr] <- rep(sum(E[, t, cr, nm])/all, all)
      
    } else if (t >= Time1) {
    
    E[Outside, t, cr, nm] <- rep(sum(E[, t - 1, cr, nm])/out, out)
    E[Inside, t, cr, nm] <- 0

    }
    
  }
  
  return(E)
  
}