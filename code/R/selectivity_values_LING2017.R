setwd("~/Projects/DensityRatio/code/R")

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'LING2017'

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

# add zero as first value in L
L1 <- append(0, L)

# manual parameters
fleet <- c('trawl', 'fixed_gear', 'WArec', 'ORrec') # names of fleets
A50up <- c(3, 5, 5, 4)
L50up <- L[A50up]
alpha <- c(0.25, 0.25, 0.55, 1)
Ffin <- c(0.07, 0, 0, 0)
A50down <- c(15, 12, 10, 9)
L50down <- L[A50down]
beta <- c(0.09, 0.3, 0.17, 0.25)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)

upcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))
downcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))

upcurve[, 1] <- downcurve[, 1] <- 0

i <- 4

upcurve[i, ages + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50up[i])))
downcurve[i, ages + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50down[i])))

plot(L1, upcurve[i, ], type = 'l', lwd = 2, col = 'black', 
     ylim = c(0, 1), 
     xlab = 'Length', ylab = 'Selectivity', 
     main = paste(species, ": ", fleet[i], sep = ''))
lines(L1, downcurve[i, ], lwd = 2, col = 'green')

# some point
abline(h = 0.5, col = 'blue')
abline(v = 56, col = 'blue')
abline(v = 87, col = 'blue')

# another point
abline(h = 0.2, col = 'purple')
abline(v = 55, col = 'purple')
abline(v = 98, col = 'purple')

# where it hits proportion = 1.0
abline(v = 59, col = 'orange')
abline(h = 1, col = 'orange')

# final value
abline(h = 0, col = 'red')