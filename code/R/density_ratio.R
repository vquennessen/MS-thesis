density_ratio <- function (t, cr, nm, A, Count, Years_sampled, Areas_sampled, 
                           Fish_sampled, Transects, Inside, Outside) {
  
  # sample all fish or just mature fish
  if (Fish_sampled == "all") {
    fish <- 1
  } else if (Fish_sampled == 'mature') {
    fish <- 2
  }
  
  # calculate count inside marine reserve
  count_in <- Count[Inside, t - 1, , fish, cr, nm]
  
  # time sampled = 1 or 3 years
  if (Years_sampled == 1) {
    years <- t - 1
  } else {
    years <- (t - Years_sampled):(t - 1)
  }
  
  # calculate counts outside marine reserve
  if (Areas_sampled == "far") {
    count_out <- Count[c(1, A), years, , fish, cr, nm]
  } else if (Areas_sampled == 'all') {
    count_out <- Count[Outside, years, , fish, cr, nm]
  }
 
  # calculate density ratio 
  DR <- (sum(count_out)/(Transects*nrow(count_out)))/(sum(count_in)/Transects)
  
  return(DR)
  
}