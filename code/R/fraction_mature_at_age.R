fraction_mature_at_age = function(Num, K_mat, L, L50) {
  
  # Initialize fraction mature at age vector
  mature <- array(rep(0, Num), c(1, Num))
  
  # Calculate fraction mature at age
  mature <- (1)/(1 + exp(K_mat*(L - L50)))
  
  return(mature)
  
}