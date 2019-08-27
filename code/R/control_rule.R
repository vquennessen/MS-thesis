control_rule <- function(a, t, cr, E, Count) {
  
  if (cr == 1) {
    E <- management(a, t, cr, E, DR = 0, CR_type = 'effort', target_DR = 0, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1)
  }
  
  if (cr == 2) {
    DR <- density_ratio(a, t, cr, Count, years_sampled = 3, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all')
    E <- management(a, t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1)
  }
  
  if (cr == 3) {
    DR <- density_ratio(a, t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all')
    E <- management(a, t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1)
  }
  
  if (cr == 4) {
    DR <- density_ratio(a, t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'far', 
                        fish_sampled = 'all')
    E <- management(a, t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1)
  }
  
  if (cr == 5) {
    DR <- density_ratio(a, t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature')
    E <- management(a, t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1)
  }
  
  if (cr == 6) {
    DR <- density_ratio(a, t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all')
    E <- management(a, t, cr, E, DR, CR_type = 'DR', target_DR = 0.8, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1)
  }
  
  if (cr == 7) {
    E <- management(a, t, cr, E, DR = 0, CR_type = 'effort', target_DR = 0, 
                    floor_DR = 0.2, effort_inc_allowed = 0)
  }
  
  if (cr == 8) {
    DR <- density_ratio(a, t, cr, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all')
    E <- management(a, t, cr, E, DR, CR_type = 'DR', target_DR = 0.8, 
                    floor_DR = 0.2, effort_inc_allowed = 0)
  }
  
  return(E)
  
}