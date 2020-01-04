################################################################################
##### Calculate alpha values
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'black rockfish 2003'

fleet <- c('sport', 'hook', 'trawl')
L50up <- c(3, 5, 10)
Sa <- c(0.4, 0.25, 0.08)
La <- c(3, 4, 8)

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
L <- length_at_age(rec_age, max_age, L1f, L2f, Kf, a1f, a2f, all_ages = F) # length at age

alpha <- (-1 * log(1/Sa - 1)) / (La - L50up)

################################################################################
##### Calculate beta values
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'black rockfish 2003'

fleet <- c('sport', 'hook', 'trawl')
Ffin <- c(0.27, 0.28, 1)
L50down <- c(7, 16, 50)
Sa <- c(0.8, 0.8, 1)
La <- c(6, 10, 50)

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
L <- length_at_age(rec_age, max_age, L1f, L2f, Kf, a1f, a2f, all_ages = F) # length at age

beta <- (-1 * log((1 - Ffin)/(1 - Sa) - 1)) / (La - L50down)

################################################################################
##### try to graph selectivity
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'black rockfish 2003'

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
L <- length_at_age(rec_age, max_age, L1f, L2f, Kf, a1f, a2f, all_ages = F) # length at age

# manual parameters
fleet <- c('sport', 'hook', 'trawl')
L50upA <- c(2, 5, 10)
L50upL <- L[L50upA - rec_age + 1]
alpha <- c(0.35, 0.6, 0.64)
Ffin <- c(0.25, 0.06, 0)
L50downA <- c(6, 16, 0)
L50downL <- L[L50downA - rec_age + 1]
beta <- c(1.2, 0.6, 0)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)

upcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))
downcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))

for (j in 1:rec_age) {
  upcurve[, j] <- downcurve[, j] <- 0
}

i <- 3

# for (i in 1:f) {
  upcurve[i, age + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50upL[i])))
  downcurve[i, age + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50downL[i])))
# }

plot(0:max_age, upcurve[i, ], type = 'l', lwd = 2, col = 'black', 
     ylim = c(0, 1))
lines(0:max_age, downcurve[i, ], lwd = 2, col = 'green')
abline(h = 0.5, col = 'blue')
abline(v = 10, col = 'blue')
abline(v = 4, col = 'orange')
abline(h = 1, col = 'red')