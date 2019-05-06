# Selectivity
# Based on Babcock & MacCall (2011): Eq. (8)
Sp1 <- (1)/(1 + exp(-1*fleet_alpha[1]*(L - L50)))
# Based on Babcock & MacCall (2011): Eq. (9)
Sp2 <- 1 - (1 - Ffin)/(1 + exp(-1*fleet_beta[1]*(L - L50)))
# TODO: Figure out what Ffin is
# TODO: Figure out selectivity, alpha and beta values

# Selectivity at Age
# Based on Babcock & MacCall (2011): Eq. (7)
Sp <- array(rep(0, n*length(fleet_alpha)), c(n, length(fleet_alpha)))