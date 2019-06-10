density_ratio <- function (a, t, count_sp, years_sampled, fished_areas_sampled, 
                        fish_sampled) {
  
  
  # error handling
  if (fish_sampled =! "all" | "mature") {
    stop("CR_type must be 'effort' or 'DR'")
  }
  
  
  if (fish_sampled == "all") {
    fish <- 1
  } else {
    fish <- 2
  }
  
  inside <- median(1:A)
  
  # calculate count inside marine reserve
  count_in <- count_sp[inside, t, , fish]
  
  # calculate counts outside marine reserve
  if (fished_areas_sampled == "far") {
    count_out <- count_sp[c(1, A), t, , fish]
  } else {count_out <- count_sp[-inside, t, , fish]}
 
  # calculate density ratio 
  DR <- (sum(count_out)/(transects*nrow(count_out)))/(sum(count_in)/transects)
  
  return(DR)
  
}