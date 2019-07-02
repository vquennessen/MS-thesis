effort_allocation <- function(a, t, allocation, A, E, biomass) {
  
  if (t < time) {
    
    if (allocation == 'IFD') {
      prop_biomass <- biomass[, t]/sum(biomass[, t])
      E <- sum(E)*prop_biomass
    }
    
  } else {
    
    areas <- 1:A
    reserve <- median(A)
    outside <- areas[-reserve]
    
    E[reserve] <- 0
    
    if (allocation == "equal") {
      E[outside] <- rep(sum(E)/(A-1), A - 1)
    } else {
      prop_biomass <- biomass[outside, t]/sum(biomass[outside, t])
      E[outside] <- sum(E)*prop_biomass
    }
    
  }
  
  return(E)
  
}