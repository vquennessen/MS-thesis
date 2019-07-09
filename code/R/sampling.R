sampling <- function(a, t, cr, r, D, abundance_all, abundance_mature, transects, 
                     x, count_sp, nuS) {
  
  # Calculate delta - constant of proportionality
  # Based on Babcock & MacCall (2011): Eq. (13)
  delta <- r / D
  
  # Total population size across all areas
  total_all <- sum(abundance_all[, t, cr])
  total_mature <- sum(abundance_mature[, t, cr])
  
  # Calculate odds ratio of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  odds_all <-  (delta * abundance_all[a, t, cr]) / (total_all / A)
  odds_mature <-  (delta * abundance_mature[a, t, cr]) * (total_mature / A)
  
  # Calculate probability based on odds ratio
  p_all <- 1 / (1 + exp(odds_all))
  p_mature <- 1 / (1 + exp(odds_mature))

  # TODO: figure out how to properly calculate probability
  
  # Determine if species is seen at least once
  # Dimensions = 1 * transects
  presence_all <- rbinom(transects, 1, p_all)
  presence_mature <- rbinom(transects, 1, p_mature)
  
  # Based on Babcock & MacCall (2011): Eq. (16)
  gamma_sp <- x / D
  
  # Calculate species count given transects with positive visuals
  count_sp[a, t, , 1, cr] <- presence_all*(gamma_sp*abundance_all[a, t, cr]*exp(nuS[a, t, cr]))
  count_sp[a, t, , 2, cr] <- presence_mature*(gamma_sp*abundance_mature[a, t, cr]*exp(nuS[a, t, cr]))
  
  return(count_sp)
}