recruitment = function(a, t, cr, nm, SSB, A, R0, h, B0, Eps, sigma_R, rec_age, 
                       recruitment_mode) {
  
  # Recruitment
  # Based on Babcock & MacCall (2011): Eq. (3)
  # Dimensions = 1 * 1
  
  if (recruitment_mode == 'closed') { ssb <- SSB[a, t - rec_age, cr, nm]}
  if (recruitment_mode == 'pool') { ssb <- sum(SSB[, t - rec_age, cr, nm] / A)}
  
  adjR0 <- R0 / A
  adjB0 <- B0 / A
  
  R1 <- (0.8 * adjR0 * h * ssb) / (0.2 * adjB0 * (1 - h) + (h - 0.2) * ssb) 
  R <- R1 * (exp(Eps[a, t, cr, nm] - sigma_R^2 / 2))

  return(R)
  
}