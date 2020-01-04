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
L <- length_at_age(rec_age, max_age, L1f, L2f, Kf, a1f, a2f, all_ages = T) # length at age
ages <- 1:max_age

# manual parameters
fleet <- c('dead', 'live', 'shore', 'ocean')
A50up <- c(4, 3, 1, 2)
L50up <- L[A50up]
alpha <- c(0.33, 0.4, 0.9, 0.35)
Ffin <- c(1, 0.7, 0.07, 1)
A50down <- c(1, 17, 3, 1)
L50down <- L[A50down]
beta <- c(0, 0.35, 0.2, 0)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)

upcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))
downcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))

upcurve[, 1] <- downcurve[, 1] <- 0

i <- 4

upcurve[i, ages + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50up[i])))
downcurve[i, ages + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50down[i])))

plot(0:max_age, upcurve[i, ], type = 'l', lwd = 2, col = 'black', 
     ylim = c(0, 1), 
     xlab = 'Age', ylab = 'Selectivity', 
     main = paste(species, ": ", fleet[i], sep = ''))
lines(0:max_age, downcurve[i, ], lwd = 2, col = 'green')

# some point
abline(h = 0.5, col = 'blue')
abline(v = 0.75, col = 'blue')
abline(v = 3, col = 'blue')

# another point
abline(h = 0.8, col = 'purple')
abline(v = 2, col = 'purple')
# abline(v = 6, col = 'purple')

# where it hits proportion = 1.0
abline(v = 1.25, col = 'orange')
abline(h = 0.825, col = 'orange')

# final value
abline(h = 0.08, col = 'red')