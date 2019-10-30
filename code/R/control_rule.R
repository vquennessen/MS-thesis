control_rule <- function(t, cr, nm, E, Count, time1, time2, transects, M) {
  
  if (cr == 1) {
    estimate <- 'low'
    DR <- density_ratio(t, cr, estimate, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    E <- management(t, cr, estimate, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 2) {
    estimate <- 'correct'
    DR <- density_ratio(t, cr, estimate, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    E <- management(t, cr, estimate, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 3) {
    estimate <- 'high'
    DR <- density_ratio(t, cr, estimate, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    E <- management(t, cr, estimate, E, DR, CR_type = 'DR', target_DR = 0.6, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 4) {
    estimate <- 'low'
    DR <- density_ratio(t, cr, estimate, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.6, 
                              estimate)
    E <- management(t, cr, nm, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 5) {
    estimate <- 'correct'
    DR <- density_ratio(t, cr, estimate, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.6, 
                              estimate)
    E <- management(t, cr, nm, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  if (cr == 6) {
    estimate <- 'correct'
    DR <- density_ratio(t, cr, estimate, Count, years_sampled = 1, 
                        fished_areas_sampled = 'all', 
                        fish_sampled = 'mature', transects)
    target_DR <- transient_DR(nat_mortality, t, time1, time2, final_DR = 0.6, 
                              estimate)
    E <- management(t, cr, nm, E, DR, CR_type = 'DR', target_DR, 
                    floor_DR = 0.2, effort_inc_allowed = 0.1, time1)
  }
  
  return(E)
  
}