setwd("~/Projects/DensityRatio/code/R")

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'CAB2019'

# load species parameters
par <- parameters(species)

max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
rec_age                <- par[[3]]        # age at recruitment
af  <- par[[4]];   bf  <- par[[5]]        # weight at length parameters (f)
a1f <- par[[6]];  L1f  <- par[[7]]        # growth parameters (f)
a2f <- par[[8]];  L2f  <- par[[9]] 
Kf  <- par[[10]]  
L50                    <- par[[11]]       # length at 50% maturity
k_mat                  <- par[[12]]       # slope of maturity curve

# Calculated values
age <- rec_age:max_age                          # applicable ages
n <- length(age)                                # number of age bins
L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f) # length at age

# manual parameters
fleet <- c('dead', 'live', 'shore', 'ocean')
A50up <- c(1, 3, 1, 2)
L50up <- L[A50up - rec_age + 1]
alpha <- c(0.4, 0.25, 0.15, 0.35)
Ffin <- c(1, 0.7, 0.43, 1)
A50down <- c(1, 17, 4, 1)
L50down <- L[A50down - rec_age + 1]
beta <- c(0, 0.35, 0.26, 0)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)

upcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))
downcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))

for (j in 1:rec_age) {
  upcurve[, j] <- downcurve[, j] <- 0
}

i <- 1

upcurve[i, age + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50up[i])))
downcurve[i, age + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50down[i])))

plot(0:max_age, upcurve[i, ], type = 'l', lwd = 2, col = 'black', 
     ylim = c(0, 1))
lines(0:max_age, downcurve[i, ], lwd = 2, col = 'green')
# some point
abline(h = 0.8, col = 'blue')
abline(v = 5, col = 'blue')
# another point
abline(h = 0.5, col = 'purple')
abline(v = 1, col = 'purple')
# abline(v = 6, col = 'purple')

# where it hits proportion = 1.0
abline(v = 12, col = 'orange')
# abline(h = 0.85, col = 'orange')
# final value
abline(h = 1, col = 'red')