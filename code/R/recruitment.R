recruitment = function(N, W, M, s, R0, h, B0, e, sigma_R) {
  
  # Spawning Stock Biomass
  B <- sum(N * W * M)
  
  # Recruitment
  # Based on Babcock & MacCall (2011): Eq. (3)
  R <- array(rep(0, s), s)
  
}