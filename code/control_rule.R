control_rule <- function(t, cr, nm, A, E, Count, Time1, TimeT, Transects, 
                         Nat_mortality, Final_DR, Inside, Outside, 
                         Areas_sampled, Fish_sampled, Years_sampled) {
  
  if (cr == 1) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 1, A, Count, Years_sampled, Areas_sampled, 
                        Fish_sampled, Transects, Inside, Outside)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = Final_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, Time1)
  }
  
  if (cr == 2) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 2, A, Count, Years_sampled, Areas_sampled,
                        Fish_sampled, Transects, Inside, Outside)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = Final_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, Time1)
  }
  
  if (cr == 3) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 3, A, Count, Years_sampled, Areas_sampled, 
                        Fish_sampled, Transects, Inside, Outside)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = Final_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, Time1)
  }
  
  target_DR <- transient_DR(start_time = Time1, end_time = timeT, Final_DR, 
                            Nat_mortality, nm)
  
  if (cr == 4) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 1, A, Count, Years_sampled, Areas_sampled, 
                        Fish_sampled, Transects, Inside, Outside)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR[t - Time1 + 1], 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, Time1)
  }
  
  if (cr == 5) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 2, A, Count, Years_sampled, Areas_sampled, 
                        Fish_sampled, Transects, Inside, Outside)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR[t - Time1 + 1], 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, Time1)
  }
  
  if (cr == 6) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 3, A, Count, Years_sampled, Areas_sampled, 
                        Fish_sampled, Transects, Inside, Outside)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR[t - Time1 + 1], 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, Time1)
  }
  
  return(E)
  
}