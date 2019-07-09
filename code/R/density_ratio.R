density_ratio <- function (a, t, cr, count_sp, years_sampled, 
                           fished_areas_sampled, fish_sampled) {
  
  # set marine reserve area as middle area
  areas <- 1:A
  inside <- median(1:A)
  y <- areas[-inside]  
  
  # sample all fish or just mature fish
  if (fish_sampled == "all") {
    fish <- 1
  } else {
    fish <- 2
  }
  
  # calculate count inside marine reserve
  count_in <- count_sp[inside, t, , fish, cr]
  
  # time sampled = 1 or 3 years
  if (years_sampled == 1) {
    years <- t
  } else {
    years <- (t - years_sampled + 1):t
  }
  
  # calculate counts outside marine reserve
  if (fished_areas_sampled == "far") {
    count_out <- count_sp[c(1, A), years, , fish, cr]
  } else {
    count_out <- count_sp[y, years, , fish, cr]
  }
 
  # calculate density ratio 
  DR <- (sum(count_out)/(transects*nrow(count_out)))/(sum(count_in)/transects)
  
  return(DR)
  
}