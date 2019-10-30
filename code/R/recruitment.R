recruitment = function(a, t, cr, SSB, A, R0, h, B0, Eps, sigma_R, rec_age) {
  
  # Recruitment
  # Based on Babcock & MacCall (2011): Eq. (3)
  # Dimensions = 1 * 1
  
  adjR0 <- R0 / A
  adjB0 <- B0 / A
  
  R1 <- (0.8 * adjR0 * h * SSB[a, t - rec_age, cr]) / (0.2 * adjB0 * (1 - h) + 
                                              (h - 0.2) * SSB[a, t - rec_age, cr]) 
  R <- R1 * (exp(Eps[a, t, cr] - sigma_R^2 / 2))

  return(R)
  
}