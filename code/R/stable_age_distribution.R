stable_age_distribution <- function(b, c, max_age, m, L0, W0, rec_age, M, Fb) {
  
  # rows for each age class, including age 0
  ages <- 0:max_age
  
  # Define egg production
  fecundity <- (b*W0 + c)*1000
  
  # Define Leslie matrix
  A <- rep(1, length(ages) - 1)
  B <- exp(-1*(M + Fb*(ages[1:length(ages) - 1] >= rec_age)))
  L1 <- diag(A*B, length(A))
  L2 <- fecundity*(ages >= m)
  
  f <- t <- matrix(rep(0, length(ages)^2), nrow = length(ages))
  t[2:nrow(t), 1:(ncol(t) - 1)] <- L1
  f[1,] <- L2
  
  LM <- t + f

  # pull out eigenvectors (V) and eigenvalues (L) of A0
  e <- eigen(LM)
  V <- Re(e$vectors)
  L0 <- e$values
  
  L <- Re(L0)
  n0 <- V[,which(L == max(L))]
  
  # dominant eigenvalue starts stable age distribution, 
  # standardized to be a unit vector
  SAD <- matrix(n0/sum(n0), ncol = 1)
  
  return(SAD)
  
}