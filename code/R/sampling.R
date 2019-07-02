sampling <- function(a, t, y, r, D, abundance_all, abundance_mature, transects, 
                     x, count_sp, nuS) {
  
  # Calculate delta - constant of proportionality
  # Based on Babcock & MacCall (2011): Eq. (13)
  delta <- r / D
  
  # Total population size across all areas
  total_all <- sum(abundance_all[, t])
  total_mature <- sum(abundance_mature[, t])
  
  # Calculate odds ratio of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  odds_all <-  (delta * abundance_all[a, t]) / (total_all / A)
  odds_mature <-  (delta * abundance_mature[a, t]) * (total_mature / A)
  
  # Calculate probability based on odds ratio
  p_all <- 1 / (1 + e^odds_all)
  p_mature <- 1 / (1 + e^odds_mature)

  # TODO: figure out how to properly calculate probability
  
  # Determine if species is seen at least once
  # Dimensions = 1 * transects
  presence_all <- rbinom(transects, 1, p_all)
  presence_mature <- rbinom(transects, 1, p_mature)
  
  # Based on Babcock & MacCall (2011): Eq. (16)
  gamma_sp <- x / D
  
  # Calculate species count given transects with positive visuals
  count_sp[a, t, , 1] <- presence_all*(gamma_sp*abundance_all[a, t]*exp(nuS[a, t, y]))
  count_sp[a, t, , 2] <- presence_mature*(gamma_sp*abundance_mature[a, t]*exp(nuS[a, t, y]))
  
  return(count_sp)
}