management <- function(a, t, E, DR, CR_type, target_DR, effort_inc_allowed) {
  
  if ((CR_type == 'DR' & DR > target_DR) | CR_type == 'effort') {
    E <- E*(1 + effort_inc_allowed)
  } else { E <- E }
  
  return(E)
  
}