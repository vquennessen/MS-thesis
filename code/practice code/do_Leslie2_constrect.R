########## 2nd Leslie function, with constant recruitment ######################

#' Run Leslie matrix model for analytical MPA stuff
#' This version: parsed from do_Leslie; initial vector n_init is specified,
#' there is no initial fished stage
#' 
#' Parsed from do_Leslie2, this version has constant recruitment
#' instead of a closed population

source("C:/Users/Vic/Box/Quennessen_Thesis/practice code/Leslie.R")

do_Leslie2_constrect = function(max_age, a_mat, t_c, L_inf, k, t0, a2, b2, a, 
                                M, t, n_init, f = 0, Rfact = 0) {
  
  # Leslie matrix, unfished
  A0 = Leslie(max_age, a_mat, t_c, L_inf, k, t0, a2, b2, a, M, 0)[[1]]
  
  # Leslie matrix, fished
  AF = Leslie(max_age, a_mat, t_c, L_inf, k, t0, a2, b2, a, M, f)[[1]]
  
  # delete reproduction
  A1 <- A0
  A1[1, ] <- 0
  
  AF1 <- AF
  AF1[1, ] <- 0
  
  # pull out eigenvectors (V) and eigenvalues (L) of A0
  e = eigen(A0)
  V = Re(e$vectors)
  L0 = e$values
  
  L = Re(L0)
  n0 = V[,which(L == max(L))]
  # dominant eigenvalue starts age distribution, standardized so n0(1) = 1
  n0 = matrix(n0/n0[1], ncol = 1)
  
  L = matrix(L0[order(Mod(L0), decreasing = T)], ncol = 1)
  L1 = L[1]
  L2 = L[2]
  
  # Now remove fishing, and run the model
  
  # initialize nmm vector
  nmm = matrix(rep(0, max_age + 1), ncol = 1)
  
  # set initial vector
  if (sum(is.nan(n_init) == T) > 0) {
    nmm = n0
  } else {
    nmm = matrix(n_init, ncol = 1)
  }
  
  R_init = rep(0, length(n0))
  R_init[1] = nmm[1]
  R_tmp = R_init
  
  for (i in seq(2, t+1)) {
    R_tmp[1] = R_tmp[1] + Rfact
    nmm = cbind(nmm, AF1 %*% nmm[, i-1] + R_tmp)
  }
  
  # measure deviation (eq. from Cohen 1979 SIAM J Appl Math)
  l1 = max(Re(eigen(A0)$values))
  
  # B = (l1 .^ -1e3) * A0 ^ 1e3
  d = eigen(A0)$vectors
  l = Re(eigen(A0)$values)
  v = matrix(abs(d[,which(l == max(l))]), ncol = 1)
  
  d = eigen(t(A0))$vectors
  l = Re(eigen(t(A0))$values)
  w = matrix(abs(d[,which(l == max(l))]), ncol = 1)
  
  B = (v%*%t(w))/as.numeric((t(v)%*%w))
  
  Z = solve(diag(ncol(B)) + B - A0/l1)
  
  D = matrix(rep(0, t + 1), ncol = t + 1)
  
  for (i in seq(1, t+1)) {
    ZB = (Z-B)%*%nmm[, i]/sum(nmm[, i])
    D[i] = sum(abs(ZB))
  }
  
  output = list(n0, nmm, L1, L2, D)
  
  return(output)
  
}
