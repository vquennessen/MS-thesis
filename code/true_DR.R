true_DR <- function(t, cr, Abundance_mature, Outside, Inside, Density_Ratios, 
                    Time1) {
  
  # Density of fish outside marine reserve(s)
  Outside_density <- sum(Abundance_mature[Outside, t, cr, 2]) / length(Outside)
  
  # Density of fish inside marine reserve(s)
  Inside_density <- sum(Abundance_mature[Inside, t, cr, 2]) / length(Inside)
  
  # True density ratio
  Density_Ratios[t - Time1 + 1, cr] <- Outside_density / Inside_density
  
  return(Density_Ratios)
}