management <- function(E, DR, CR_type, target_DR, floor_DR, effort_inc_allowed) {
  
  if (CR_type == 'effort') {
    
    E <- E*(1 + effort_inc_allowed)
    
  }
  
  if (DR > target_DR) {
    
    E <- E*(1 + effort_inc_allowed)
    
  }
  
  if (DR < target_DR && DR > floor_DR) {
    
    E <- E*(1 - effort_inc_allowed)
    
  }
  
  if (DR < floor_DR) {
    
    E <- E*0
    
  }
  
  if (sum(E) == 0) {
    E <- rep(0.01, A)
  }
  
  return(E)
  
}