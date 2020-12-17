# clear environment
rm(list = ls())

# source required functions
remotes::install_github('vquennessen/densityratio')
library(densityratio)

species <- 'China2015'

# # load species parameters
# par <- parameters(species)

max_age                 <- 79       # maximum age
M                       <- 0.07     # natural mortality
rec_age                 <- 5        # age at recruitment
a  <- 7.79E-06;   b     <- 3.177    # weight at length parameters (f)
A1 <- 10;         L1    <- 30.5        # growth parameters (f)
A2 <- 30;         L2    <- 36.85 
K  <- 0.147  
L50                     <- 27       # length at 50% maturity
k_mat                   <- -0.467   # slope of maturity curve

# Calculated values
L <- length_age(rec_age, max_age, A1, L1, A2, L2, K, All_ages = T) # length at age
ages <- 1:max_age

# add zero as first value in L
L1 <- append(0, L)

# manual parameters
fleet <- c('Commercial', 'Rec_PC', 'Rec_PR') # names of fleets
A50up <- c(10, 5, 5)
L50up <- L[A50up]
alpha <- c(1, 0.25, 0.55)
Ffin <- c(1, 1, 1)
A50down <- c(79, 79, 79)
L50down <- L[A50down]
beta <- c(1, 1, 1)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)

upcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))
downcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))

upcurve[, 1] <- downcurve[, 1] <- 0

i <- 1

upcurve[i, ages + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50up[i])))
downcurve[i, ages + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50down[i])))

plot(L1, upcurve[i, ], type = 'l', lwd = 2, col = 'black', 
     ylim = c(0, 1), 
     xlab = 'Length', ylab = 'Selectivity', 
     main = paste(species, ": ", fleet[i], sep = ''))
lines(L1, downcurve[i, ], lwd = 2, col = 'green')

# some point before the concavity switches
abline(h = 0.275, col = 'blue')
abline(v = 30, col = 'blue')

# another point after the concavity switches
abline(h = 0.825, col = 'purple')
abline(v = 32.5, col = 'purple')

# where it hits proportion = 1.0
abline(v = 35, col = 'orange')
abline(h = 1, col = 'orange')
