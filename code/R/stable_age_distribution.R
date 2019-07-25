stable_age_distribution <- function(b, c, max_age, m, L0, W0, rec_age, M, Fb) {
  
  # rows for each age class, including age 0
  ages <- 0:max_age
  
  # Define egg production
  fecundity <- (b*W0 + c)*1000
  
  # Define Leslie matrix
  Y <- rep(1, max_age)
  Z <- exp(-1*(M + Fb*(ages[1:max_age] >= rec_age)))
  L1 <- diag(Y*Z, length(Y))
  L2 <- fecundity*(ages >= m)
  
  one <- two <- matrix(rep(0, length(ages)^2), nrow = length(ages))
  one[2:nrow(one), 1:(ncol(one) - 1)] <- L1
  two[1, ] <- L2
  
  LM <- one + two

  # pull out eigenvectors (V) and eigenvalues (L) of A0
  e <- eigen(LM)
  V <- Re(e$vectors)
  L <- Re(e$values)
  
  n0 <- V[,which(L == max(L))]
  
  # dominant eigenvalue starts stable age distribution, 
  # take out ages before recruitment
  SAD <- n0[(rec_age + 1):(max_age + 1)]
  
  return(SAD)
  
}
