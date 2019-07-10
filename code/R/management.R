management <- function(a, t, cr, E, DR, CR_type, target_DR, floor_DR, effort_inc_allowed) {
  
  if (sum(E[, t - 1, cr]) == 0) {
    
    E[, t, cr] <- E[, time, cr]*0.10
    
  } else if (CR_type == 'effort' | DR > target_DR) {
    
    E[a, t, cr] <- E[a, t - 1, cr]*(1 + effort_inc_allowed)
    
  } else if (DR < target_DR & DR > floor_DR) {
    
    E[a, t, cr] <- E[a, t - 1, cr]*(1 - effort_inc_allowed)
    
  } else if (DR < floor_DR && sum(E[a, t, cr]) != 0) {
    
    E[a, t, cr] <- E[a, t - 1, cr]*0
    
  }
  
  return(E)
  
}