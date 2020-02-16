##### start #####

setwd("~/Projects/MS-thesis/code")

# clear environment
rm(list = ls())

# source required package
devtools::install_git('https://github.com/vquennessen/densityratio.git', 
                      method = 'libcurl', force = TRUE)
library(densityratio)

species <- 'ES_COW_2013'

# load species parameters
par <- parameters(Species = 'ES_COW_2013')

Max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
Rec_age                <- par[[3]]        # age at recruitment
WA  <- par[[4]];  WB   <- par[[5]]        # weight at length parameters (f)
A1  <- par[[6]];  L1   <- par[[7]]        # growth parameters (f)
A2  <- par[[8]];  L2   <- par[[9]] 
K  <- par[[10]]  
L50                    <- par[[11]]       # length at 50% maturity
# k_mat                  <- par[[12]]       # slope of maturity curve

# Length at age
L <- length_age(Rec_age, Max_age, A1, L1, A2, L2, K, All_ages = T)
ages <- 1:Max_age
plot(ages, L)
# add zero as first value in L
L1 <- append(0, L)

# Maturity at age (estimate K_mat)
K_mat <- -0.25
Mat <- maturity(Rec_age, Max_age, K_mat, L, L50)
plot(L[1:40], Mat[1:40], ylim = c(0, 1))
points(x = c(38, 42.5, 47.5), y = c(0.2, 0.5, 0.8), pch = 4, col = 'blue')

##### manual parameters ####
fleet <- c('trawl', 'non_trawl', 'rec', 'hake', 'research') # names of fleets
A50up <- c(5, 5, 4, 8, 1)
L50up <- L[A50up]
alpha <- c(0.3, 0.6, 1, 1, 1)
Ffin <- c(0.36, 1, 0.175, 0.65, 0.8)
A50down <- c(10, 50, 7, 11, 30)
L50down <- L[A50down]
beta <- c(1, 0, 1, 1, 0.08)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)

upcurve <- array(rep(NA, f * (Max_age + 1)), c(f, Max_age + 1))
downcurve <- array(rep(NA, f * (Max_age + 1)), c(f, Max_age + 1))

upcurve[, 1] <- downcurve[, 1] <- 0

i <- 5

upcurve[i, ages + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50up[i])))
downcurve[i, ages + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50down[i])))

plot(L1, upcurve[i, ], type = 'l', lwd = 2, col = 'black', 
     ylim = c(0, 1), 
     xlab = 'Length', ylab = 'Selectivity', 
     main = paste(species, ": ", fleet[i], sep = ''))
lines(L1, downcurve[i, ], lwd = 2, col = 'green')

##### values from ss file #####
lsel2 <- c(5.91436e-005, 8.83504e-005, 0.000123398, 0.00252993, 0.00718828,
           0.0157125, 0.030458, 0.0545767, 0.0919005, 0.146575, 0.222376, 
           0.321622, 0.443617, 0.582771, 0.727095, 0.858337, 0.955035,
           0.998426, 0.999979, 0.806091, 0.428221, 0.362598, 0.360571,
           0.360559, 0.360559, 0.360559, 0.360559, 0.360559, 0.360559, 0.360559)

lsel5 <- c(5.91438e-005, 8.83506e-005, 0.000123399, 0.000241094, 0.000663171,
           0.00201613, 0.00588791, 0.015764, 0.038171, 0.0832551, 0.16336,
           0.288235, 0.457243, 0.652104, 0.836076, 0.963668, 0.99999, 1, 0.999868,
           0.999452, 0.998763, 0.997823, 0.996661, 0.995312, 0.993812, 0.992204,
           0.990528, 0.988823, 0.987127, 0.985474)

lsel11 <- c(5.91436e-005, 8.83503e-005, 0.000123398, 0.000123398, 0.000123399,
            0.000123399, 0.000123399, 0.000123399, 0.0001234, 0.000123416,
            0.000123974, 0.000137379, 0.000351778, 0.00262228, 0.0184335,
            0.0899669, 0.295336, 0.649703, 0.957279, 0.999994, 0.853148,
            0.647117, 0.646909, 0.646909, 0.646909, 0.646909, 0.646909,
            0.646909, 0.646909, 0.646909)

lsel17 <- c(0.199126, 0.29746, 0.41546, 0.978568, 0.994456, 0.996059, 0.996005,
            0.995311, 0.994188, 0.992695, 0.990847, 0.988651, 0.986117, 0.983252,
            0.980069, 0.976582, 0.972806, 0.968753, 0.96444, 0.959879, 0.955087,
            0.950077, 0.944865, 0.939466, 0.933895, 0.928168, 0.9223, 0.916306,
            0.910204, 0.904007)

##### plotting aids ######

# points
points(x = seq(9, 67, by = 2), y = lsel17)

# # some point
# abline(h = 0.5, col = 'blue')
# abline(v = 27, col = 'blue')
# abline(v = 40, col = 'blue')
# 
# # another point
# abline(h = 0.8, col = 'purple')
# abline(v = 29.5, col = 'purple')
# abline(v = 37, col = 'purple')
# 
# # where it hits proportion = 1.0
# abline(v = 34, col = 'orange')
# abline(h = 1, col = 'orange')
# 
# # final value
# abline(h = 0.175, col = 'red')