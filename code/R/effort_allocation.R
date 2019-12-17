effort_allocation <- function(t, cr, nm, allocation, E, yield, time1, inside, 
                              outside) {
  
  #' number of areas not in a reserve
  out <- length(outside)
  ins <- length(inside)
  all <- out + ins
  
  #' If effort is allocated using the ideal free distribution, effort for one 
  #' year depends on the distribution of yield from the previous year
  if (allocation == 'IFD') {
      
      prop_yield_outside <- yield[outside, t - 1, cr, nm]/sum(yield[outside, t - 1, cr, nm])
      
      E[outside, t, cr, nm] <- sum(E[, t - 1, cr, nm])*prop_yield_outside
      E[inside, t, cr, nm] <- 0
    
  #' Otherwise, distribute effort equally between the four areas outside the 
  #' marine reserve, regardless of yield
  } else if (allocation == 'equal') {
    
    if (t < time1) {
      
      E[, t, cr] <- rep(sum(E[, t, cr, nm])/all, all)
      
    } else if (t >= time1) {
    
    E[outside, t, cr, nm] <- rep(sum(E[, t - 1, cr, nm])/out, out)
    E[inside, t, cr, nm] <- 0

    }
    
  }
  
  return(E)
  
}