effort_allocation <- function(t, cr, nm, allocation, A, E, yield, time1) {
  
  # separate areas inside and outside of reserve
  areas <- 1:A
  reserve <- median(areas)
  outside <- areas[-reserve]
  
  #' If effort is allocated using the ideal free distribution, effort for one 
  #' year depends on the distribution of yield from the previous year
  if (allocation == 'IFD') {
      
      prop_yield_outside <- yield[outside, t - 1, cr, nm]/sum(yield[outside, t - 1, cr, nm])
      
      E[outside, t, cr, nm] <- sum(E[, t - 1, cr, nm])*prop_yield_outside
      E[reserve, t, cr, nm] <- 0
    
  #' Otherwise, distribute effort equally between the four areas outside the 
  #' marine reserve, regardless of yield
  } else if (allocation == 'equal') {
    
    if (t < time1) {
      
      E[, t, cr] <- rep(sum(E[, t, cr, nm])/A, A)
      
    } else if (t >= time1) {
    
    E[reserve, t, cr, nm] <- 0
    E[outside, t, cr, nm] <- rep(sum(E[, t - 1, cr, nm])/(A - 1), A - 1)
    
    }
    
  }
  
  return(E)
  
}