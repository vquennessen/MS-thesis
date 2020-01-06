fraction_mature_at_age = function(n, k_mat, L, L50) {
  
  M <- array(rep(0, n), c(1, n))
  M <- (1)/(1 + exp(k_mat*(L - L50)))
  
  return(M)
  
}