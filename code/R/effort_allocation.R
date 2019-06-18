effort_allocation <- function(allocation, E, biomass, a, t) {
  
  if (allocation == "equal") {
    E <- rep(sum(E)/A, A)
  } else {
    prop_biomass <- biomass[a, t]/sum(biomass[a, t])
    E <- sum(E)*prop_biomass
  }
  
}