effort_allocation <- function(t, cr, allocation, A, E, yield, time1) {
  
  # separate areas inside and outside of reserve
  areas <- 1:A
  reserve <- median(areas)
  outside <- areas[-reserve]
  
  #' If effort is allocated using the ideal free distribution, effort for one 
  #' year depends on the distribution of yield from the previous year
  if (allocation == 'IFD') {
    
    # Before reserve implementation
    if (t < time1) {
      
      prop_yield <- yield[, t - 1, cr]/sum(yield[, t - 1, cr])
      E[, t, cr] <- sum(E[, t - 1, cr])*prop_yield
      
      #' After reserve implementation, effort in the reserve is equal to zero
    } else if (t >= time1) {
      
      prop_yield_outside <- yield[outside, t - 1, cr]/sum(yield[outside, t - 1, cr])
      
      E[outside, t, cr] <- sum(E[, t - 1, cr])*prop_yield_outside
      E[reserve, t, cr] <- 0
      
    }
    
  #' Otherwise, distribute effort equally between the four areas outside the 
  #' marine reserve, regardless of yield
  } else if (allocation == 'equal') {
    
    if (t < time1) {
      
      E[, t, cr] <- rep(sum(E[, t, cr])/A, A)
      
    } else if (t >= time1) {
    
    E[reserve, t, cr] <- 0
    E[outside, t, cr] <- rep(sum(E[, t-1, cr])/(A-1), A-1)
    
    }
    
  }
  
  return(E)
  
}