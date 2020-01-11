sampling <- function(a, t, cr, nm, Delta, Gamma, Abundance_all, 
                     Abundance_mature, Transects, X, Count, NuS, A) {
  
  # Total population size across all areas
  total_all <- sum(Abundance_all[, t, cr, nm])
  total_mature <- sum(Abundance_mature[, t, cr, nm])
  
  # Calculate odds ratio of seeing a fish
  # Based on Babcock & MacCall (2011): Eq. (12)
  odds_all <-  (Delta * Abundance_all[, t, cr, nm]) / (total_all / A)
  odds_mature <-  (Delta * Abundance_mature[, t, cr, nm]) / (total_mature / A)
  
  # Calculate probability based on odds ratio
  p_all <- 1 / (1 + exp(odds_all))
  p_mature <- 1 / (1 + exp(odds_mature))

  # Determine if species is seen at least once
  # Dimensions = 1 * transects
  presence_all <- rbinom(Transects, 1, p_all)
  presence_mature <- rbinom(Transects, 1, p_mature)
  
  # Calculate species count given transects with positive visuals
  Count[a, t, , 1, cr, nm] <- presence_all*(Gamma*Abundance_all[a, t, cr, nm]*exp(NuS[a, t, cr, nm]))
  Count[a, t, , 2, cr, nm] <- presence_mature*(Gamma*Abundance_mature[a, t, cr, nm]*exp(NuS[a, t, cr, nm]))
  
  return(Count[a, t, , , cr, nm])
}
