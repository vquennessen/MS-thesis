##################### Leslie Matrix Function ###################################

Leslie = function(max_age, a_mat, t_c, L_inf, k, t0, a2, b2, a, M, f) {
  
  m_at_age = 0:max_age
  
  # Define length 
  l_a = exp(-1*M*m_at_age)
  len_a = L_inf*(1-exp(-1*k*(m_at_age - t0)))
  len_a = ifelse(len_a < 0, 0, len_a) # replaces any negative len_a values with 0

  # Define size
  size_a = a2*len_a^b2
  
  # Define egg production
  egg_a = a*size_a*(m_at_age >= a_mat)
  
  
  # Define leslie matrix
  A = rep(1, length(m_at_age) - 1)
  B = exp(-1*(M + f*(m_at_age[1:length(m_at_age) - 1] >= t_c)))
  L1 = diag(A*B, length(A))
  L2 = egg_a*(m_at_age >= a_mat)
  
  F = T = matrix(rep(0, length(m_at_age)^2), nrow = length(m_at_age))
  T[2:nrow(T), 1:(ncol(T) - 1)] = L1
  F[1,] = L2
  
  L = T + F
  
  # SD of spawning age distribution
  x = l_a*egg_a*(m_at_age >= a_mat)
  s_a = m_at_age*x / sum(x)
  
  mean_age = sum(s_a)
  st_a = (m_at_age - mean_age)^2*(s_a/m_at_age)
  st_age = sqrt(sum(st_a, na.rm = TRUE))/mean_age
  
  output = list(L, mean_age, st_age)

  return(output)
}

# parameters
max_age = 35
a_mat = 8
t_c = 2
L_inf = 35
k = 0.2022
t0 = 1
a2 = 15
b2 = 3
a = 1.68*10^(-5)
M = 0.14
f = 0.2

Leslie(max_age, a_mat, t_c, L_inf, k, t0, a2, b2, a, M, f)
  
# matlab parameters
maxage = 35
amat = 8
tc = 2
Linf = 35
k = 0.2022
tzero = 1
a2 = 15
b2 = 3
a = 1.68*10^(-5)
M = 0.14
f = 0.2