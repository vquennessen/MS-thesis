control_rule <- function(a, t, E, count_sp, x) {
  
  if (x == 1) {
    E <- management(a, t, E, DR = NULL, CR_type = "effort", target_DR = NULL, 
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 2) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 3, 
                        fished_areas_sampled = "all", 
                        fish_sampled = "all")
    E <- management(a, t, E, DR, CR_type = "DR", target_DR = 0.6, 
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 3) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = "all", 
                        fish_sampled = "all")
    E <- management(a, t, E, DR, CR_type = "DR", target_DR = 0.6, 
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 4) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = "far", 
                        fish_sampled = "all")
    E <- management(a, t, E, DR, CR_type = "DR", target_DR = 0.6, 
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 5) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = "all", 
                        fish_sampled = "mature")
    E <- management(a, t, E, DR, CR_type = "DR", target_DR = 0.6, 
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 6) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = "all", 
                        fish_sampled = "all")
    E <- management(a, t, E, DR, CR_type = "DR", target_DR = 0.8, 
                    effort_inc_allowed = 0.1)
  }
  
  if (x == 7) {
    E <- management(a, t, E, DR = NULL, CR_type = "effort", target_DR = NULL, 
                    effort_inc_allowed = 0)
  }
  
  if (x == 8) {
    DR <- density_ratio(a, t, count_sp, years_sampled = 1, 
                        fished_areas_sampled = "all", 
                        fish_sampled = "all")
    E <- management(a, t, E, DR, CR_type = "DR", target_DR = 0.8, 
                    effort_inc_allowed = 0)
  }
  
  return(E)
  
}