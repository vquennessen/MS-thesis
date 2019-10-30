control_rule <- function(t, cr, E, Count, time1, time2, transects, M) {
  
  if (cr == 1) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 2) {
    E <- management(t, cr, E, DR = 0, CR_type = 'effort', target_DR = NA, 
                    floor_DR = NA, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 3) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(M, t, time1, time2, final_DR = 0.6, 
                              estimate = 'low')
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 4) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(M, t, time1, time2, final_DR = 0.6, 
                              estimate = 'correct')
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 5) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(M, t, time1, time2, final_DR = 0.6, 
                              estimate = 'high')
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  return(E)
  
}