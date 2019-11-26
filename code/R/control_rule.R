control_rule <- function(t, cr, nm, A, E, Count, time1, time2, transects, 
                         nat_mortality, Density_Ratios) {
  
  if (cr == 1) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 1, A, Count,
                        years_sampled = 1, fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.8, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 2) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 2, A, Count, 
                        years_sampled = 1, fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.8, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 3) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 3, A, Count,
                        years_sampled = 1, fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    
    # calculate effort
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.8, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 4) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 1, A, Count,
                        years_sampled = 1, fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    
    # calculate effort
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.8, nm)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 5) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 2, A, Count,
                        years_sampled = 1, fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    
    # calculate effort
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.8, nm)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 6) {
    # calculate density ratio
    DR <- density_ratio(t, cr, nm = 3, A, Count,
                        years_sampled = 1, fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    
    # calculate effort
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.8, nm)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  output <- list(E, Density_Ratios)
  
  return(output)
  
}