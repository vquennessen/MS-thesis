management <- function(t, cr, estimate, E, DR, CR_type, target_DR, floor_DR, 
                       effort_inc_allowed, time1) {
  
  nm <- ifelse(estimate == 'low', 1, ifelse(estimate == 'correct', 2, 3))
  
  #' If the control rule is based on effort, and the density ratio is higher
  #' than the target density ratio, allow effort in each area to increase by the
  #' allowed effort increase value (typically 10%) 
  if (CR_type == 'effort' || DR > target_DR) {
    
    E[, t, cr, 2] <- E[, t - 1, cr, 2]*(1 + effort_inc_allowed)
   
  #' otherwise, if the density ratio is lower than the target density ratio but
  #' greater than the floor density ratio, allow effort in each area to decrease
  #' by the allowed effort increase value (typically 10%)   
  } else if (DR < target_DR & DR > floor_DR) {
    
    E[, t, cr, 2] <- E[, t - 1, cr, 2]*(1 - effort_inc_allowed)
    
  #' Finally, if the density ratio is below the floor density ratio, effort is 
  #' decreased back down to 10% of the original value
  } else if (DR < floor_DR) {
    
    E[, t, cr, 2] <- E[, time1, cr, 2]*0.10
    
  }
  
  return(E)
  
}