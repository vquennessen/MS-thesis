control_rule <- function(t, cr, E, Count, time1, transects) {
  
  if (cr == 1) {
    E <- management(t, cr, E, DR = 0, CR_type = 'effort', target_DR = 0, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 2) {
    DR <- density_ratio(t, cr, Count, years_sampled = 3, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 3) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 4) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'far', 
                        fish_sampled = 'all', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 5) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 6) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.8, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 7) {
    E <- management(t, cr, E, DR = 0, CR_type = 'effort', target_DR = 0, 
                    floor_DR = 0.2, effort_inc_allowed = 0, time1)
  }
  
  if (cr == 8) {
    DR <- density_ratio(t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.8, 
                    floor_DR = 0.2, effort_inc_allowed = 0, time1)
  }
  
  return(E)
  
}