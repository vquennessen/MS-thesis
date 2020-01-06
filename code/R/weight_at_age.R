weight_at_age = function(L, af, bf) {
  
  # Weight at age
  # Based on Babcock & MacCall (2011): Eq. (11)
  W <- af*L^bf
  
  return(W)
  
}