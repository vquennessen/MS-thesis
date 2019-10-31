control_rule <- function(t, cr, nm, E, Count, time1, time2, transects, 
                         nat_mortality) {
  
  if (cr == 1) {
    DR <- density_ratio(t, cr, nm = 1, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 2) {
    DR <- density_ratio(t, cr, nm = 2, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 3) {
    DR <- density_ratio(t, cr, nm = 3, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 4) {
    DR <- density_ratio(t, cr, nm = 1, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.6, nm)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 5) {
    DR <- density_ratio(t, cr, nm = 2, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.6, nm)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 6) {
    DR <- density_ratio(t, cr, nm = 3, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.6, nm)
    E <- management(t, cr, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  return(E)
  
}