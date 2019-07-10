effort_allocation <- function(a, t, cr, allocation, A, E, biomass) {
  
  if (t < time) {
    
    if (allocation == 'IFD') {
      prop_biomass <- biomass[, t]/sum(biomass[, t])
      E[a, t, cr] <- sum(E[, t, cr])*prop_biomass
    }
    
  } else {
    
    areas <- 1:A
    reserve <- median(A)
    outside <- areas[-reserve]
    
    E[reserve, t, cr] <- 0
    
    if (allocation == "equal") {
      E[outside, t, cr] <- rep(sum(E[, t, cr])/(A-1), A - 1)
    } else {
      prop_biomass <- biomass[outside, t]/sum(biomass[outside, t])
      E[outside, t, cr] <- sum(E[, t, cr])*prop_biomass
    }
    
  }
  
  return(E)
  
}