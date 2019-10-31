sampling <- function(a, t, cr, nm, Delta, Gamma, abundance_all, 
                     abundance_mature, transects, x, Count, nuS) {
  
  # Total population size across all areas
  total_all <- sum(abundance_all[, t, cr, nm])
  total_mature <- sum(abundance_mature[, t, cr, nm])
  
  # Calculate odds ratio of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  odds_all <-  (Delta * abundance_all[, t, cr, nm]) / (total_all / A)
  odds_mature <-  (Delta * abundance_mature[, t, cr, nm]) / (total_mature / A)
  
  # Calculate probability based on odds ratio
  p_all <- 1 / (1 + exp(odds_all))
  p_mature <- 1 / (1 + exp(odds_mature))

  # Determine if species is seen at least once
  # Dimensions = 1 * transects
  presence_all <- rbinom(transects, 1, p_all)
  presence_mature <- rbinom(transects, 1, p_mature)
  
  # Calculate species count given transects with positive visuals
  Count[a, t, , 1, cr, nm] <- presence_all*(Gamma*abundance_all[a, t, cr, nm]*exp(nuS[a, t, cr, nm]))
  Count[a, t, , 2, cr, nm] <- presence_mature*(Gamma*abundance_mature[a, t, cr, nm]*exp(nuS[a, t, cr, nm]))
  
  return(Count[a, t, , , cr, nm])
}
