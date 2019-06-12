stable_age_distribution <- function(max_age, m, L0, W0, rec_age, af, M, Fb) {
  
    ages = 0:max_age

    # Define egg production
    egg_a = af*W0*(ages >= m)
    
    # Define leslie matrix
    A = rep(1, length(ages) - 1)
    B = exp(-1*(M + Fb*(ages[1:length(ages) - 1] >= rec_age)))
    L1 = diag(A*B, length(A))
    L2 = egg_a*(ages >= m)
    
    f = t = matrix(rep(0, length(ages)^2), nrow = length(ages))
    t[2:nrow(t), 1:(ncol(t) - 1)] = L1
    f[1,] = L2
    
    L = t + f
    
    SAD <- eigen(L)
    
    output = list(L, mean_age, st_age)
    
    return(output)
  
}