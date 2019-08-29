effort_allocation <- function(a, t, cr, allocation, A, E, biomass, time1) {
  
  if (allocation == 'IFD') {
    
    prop_biomass <- biomass[, t - 1]/sum(biomass[, t - 1, cr])
    E[a, t, cr] <- sum(E[, t - 1, cr])*prop_biomass
    
  } else if (t >= time1) {
    
    areas <- 1:A
    reserve <- median(areas)
    outside <- areas[-reserve]
    
    E[reserve, t, cr] <- 0
    
    E[outside, t, cr] <- rep(sum(E[, t - 1, cr])/(A-1), A - 1)
    
  }
  
  return(E)
  
}