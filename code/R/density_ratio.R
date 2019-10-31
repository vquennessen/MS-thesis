density_ratio <- function (t, cr, nm, Count, years_sampled, 
                           fished_areas_sampled, fish_sampled, transects) {
  
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
  count_in <- Count[inside, t, , fish, cr, nm]
  
  # time sampled = 1 or 3 years
  if (years_sampled == 1) {
    years <- t
  } else {
    years <- (t - years_sampled + 1):t
  }
  
  # calculate counts outside marine reserve
  if (fished_areas_sampled == "far") {
    count_out <- Count[c(1, A), years, , fish, cr, nm]
  } else {
    count_out <- Count[y, years, , fish, cr, nm]
  }
 
  # calculate density ratio 
  DR <- (sum(count_out)/(transects*nrow(count_out)))/(sum(count_in)/transects)
  
  return(DR)
  
}