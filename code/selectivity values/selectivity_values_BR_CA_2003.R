# clear environment
rm(list = ls())

species <- 'BR_CA_2003'

# load species parameters
par <- parameters(species)

Max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
Rec_age                <- par[[3]]        # age at recruitment
WA  <- par[[4]];   WB  <- par[[5]]        # weight at length parameters (f)
A1  <- par[[6]];   L1  <- par[[7]]        # growth parameters (f)
A2  <- par[[8]];   L2  <- par[[9]] 
K   <- par[[10]]  
L50                    <- par[[11]]       # length at 50% maturity
K_mat                  <- par[[12]]       # slope of maturity curve

# Calculated values
age <- Rec_age:Max_age                          # applicable ages
n <- length(age)                                # number of age bins
L <- length_age(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = T)
# add zero as first value in L
L1 <- append(0, L)

# manual parameters
fleet <- c('sport', 'hook', 'trawl')
L50upA <- c(2, 5, 10)
L50upL <- L[L50upA]
alpha <- c(0.35, 0.6, 0.64)
Ffin <- c(0.25, 0.06, 1)
L50downA <- c(6, 16, Max_age)
L50downL <- L[L50downA]
beta <- c(1.2, 0.6, 0)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)

upcurve <- array(rep(NA, f * (Max_age + 1)), c(f, Max_age + 1))
downcurve <- array(rep(NA, f * (Max_age + 1)), c(f, Max_age + 1))

for (j in 1:Rec_age) {
  upcurve[, j] <- downcurve[, j] <- 0
}

i <- 3

upcurve[i, age + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L[Rec_age:Max_age] - L50upL[i])))
downcurve[i, age + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L[Rec_age:Max_age] - L50downL[i])))

plot(L1, upcurve[i, ], type = 'l', lwd = 2, col = 'black', 
     ylim = c(0, 1), xlab = 'Length (cm)', ylab = 'Selectivity', 
     main = paste(species, ": ", fleet[i], sep = ""))
lines(L1, downcurve[i, ], lwd = 2, col = 'green')
abline(h = 0.5, col = 'blue')
abline(v = 10, col = 'blue')
abline(v = 4, col = 'orange')
abline(h = 1, col = 'red')
