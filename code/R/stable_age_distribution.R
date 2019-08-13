stable_age_distribution <- function(b, c, max_age, m, L0, W0, rec_age, M, Fb, 
                                    h, R0, W) {
  
  # rows for each age class, including age 0
  ages <- 0:max_age
  
  # Define egg production
  fecundity <- (b*W0 + c)*1000
  
  # Define Leslie matrix
  Y <- rep(1, max_age)
  Z <- exp(-1*(M + Fb*(ages[1:max_age] >= rec_age)))
  L1 <- diag(Y*Z, length(Y))
  L2 <- fecundity*(ages >= m)
  
  ### Multiply fecundity by "a" value, as outlined in Hilborn paper

  # B is the inverse of the asymptotic recruitment in the Beverton-Holt
  # stock recruitment curve
  B <- (h - 0.2) / (0.8 * h * R0)

  # Q1 as outlined in Lawson & Hilborn 1985
  sum1 <- 0
  for (i in rec_age:max_age) {
    sum1 <- sum1 + W[i - rec_age + 1]*exp(-M * (i - rec_age))
  }

  Q1 <- W[1] + sum1

  E0 <- R0 * Q1

  # A is the inverse of the initial slope of the Beverton-Holt stock
  # recruitment curve
  A <- Q1 - B*E0

  L2 <- L2*A
  
  ### Put matrix together
  one <- two <- matrix(rep(0, length(ages)^2), nrow = length(ages))
  one[2:nrow(one), 1:(ncol(one) - 1)] <- L1
  two[1, ] <- L2
  
  LM <- one + two

  ### pull out eigenvectors (V) and eigenvalues (L) of A0
  e <- eigen(LM)
  V <- Re(e$vectors)
  L <- Re(e$values)
  
  n0 <- V[, which(L == max(L))]
  
  # dominant eigenvalue starts stable age distribution, 
  # take out ages before recruitment
  SAD <- abs(n0[(rec_age + 1):(max_age + 1)])
  
  return(SAD)
  
}
