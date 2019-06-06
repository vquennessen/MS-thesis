management <- function (a, t, count_sp, type = "DR", target_DR, 
                        effort_inc_allowed, years_sampled, fished_areas_sampled, 
                        fish_sampled) {
  
  
  # error handling
  if (CR_type =! "effort" | "DR") {
    stop("CR_type must be 'effort' or 'DR'")
  }
  
  if (fish_sampled =! "all" | "mature") {
    stop("CR_type must be 'effort' or 'DR'")
  }
  
  if (fish_sampled == "all") {
    fish <- 1
  } else {
    fish <- 2
  }
  
  # calculate density ratio
  count_in <- count_sp[3, t, , fish]
  
  
  
  if (CR_type == "effort") {
    E <- E*effort_inc_allowed
  } else {
    
  }
  
}