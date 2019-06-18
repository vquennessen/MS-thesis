control_rule <- function(a, t, E, count_sp, x) {
  
  if (x == 1) {
    E <- management(E, DR = 0, CR_type = 'effort', target_DR = 0, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1)
  }
  
  if (x == 2) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 3, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all')
    E <- management(E, DR, CR_type = 'DR', target_DR = 0.6, floor_DR = 0.2,
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 3) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all')
    E <- management(E, DR, CR_type = 'DR', target_DR = 0.6, floor_DR = 0.2,
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 4) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = 'far', 
                        fish_sampled = 'all')
    E <- management(E, DR, CR_type = 'DR', target_DR = 0.6, floor_DR = 0.2,
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 5) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature')
    E <- management(E, DR, CR_type = 'DR', target_DR = 0.6, floor_DR = 0.2,
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 6) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all')
    E <- management(E, DR, CR_type = 'DR', target_DR = 0.8, floor_DR = 0.2,
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 7) {
    E <- management(E, DR = 0, CR_type = 'effort', target_DR = 0, 
                    floor_DR = 0.2, effort_inc_allowed = 0)
  }
  
  if (x == 8) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'all')
    E <- management(E, DR, CR_type = 'DR', target_DR = 0.8, floor_DR = 0.2,
                    effort_inc_allowed = 0)
  }
  
  return(E)
  
}