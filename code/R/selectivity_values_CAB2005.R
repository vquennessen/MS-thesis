################################################################################
##### Calculate alpha values
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'CAB2005'

fleet <- c('dead', 'live', 'man.made', 'shore', 'PBR', 'CPFV')
A50up <- c(2, 3, 2, 2, 1)
Sa <- c(0.1, 0.1, 0.8, 0.1, 0.85)
Aa <- c(1, 2, 3, 1, 2)

# all 0 values for man-made fleet

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

alpha <- (-1 * log(1/Sa - 1)) / (L[Aa] - L[A50up])

################################################################################
##### Calculate beta values
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'CAB2005'

fleet <- c('dead', 'live', 'man.made', 'shore', 'PBR', 'CPFV')
Ffin <- c(0.75, 0.38, 0.43)
A50down <- c(17, 6, 6)
Sa <- c(0.9, 0.8, 0.6)
Aa <- c(15, 4, 5)

# all values 0 for dead, PBR, and CPFV

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

beta <- (-1 * log((1 - Ffin)/(1 - Sa) - 1)) / (L[Aa] - L[A50down])

################################################################################
##### try to graph selectivity
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'CAB2005'

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
fleet <- c('dead', 'live', 'man.made', 'shore', 'PBR', 'CPFV')
A50up <- c(2, 3, 1, 1, 2, 1)
L50up <- L[A50up - rec_age + 1]
alpha <- c(0.4, 0.25, 5, 0.15, 0.35, 0.25)
Ffin <- c(1, 0.7, 0.38, 0.43, 1, 1)
A50down <- c(1, 17, 5, 4, 1, 1)
L50down <- L[A50down - rec_age + 1]
beta <- c(0, 0.35, 0.5, 0.26, 0, 0)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)

upcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))
downcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))

for (j in 1:rec_age) {
  upcurve[, j] <- downcurve[, j] <- 0
}

i <- 6

upcurve[i, age + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50up[i])))
downcurve[i, age + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50down[i])))

plot(0:max_age, upcurve[i, ], type = 'l', lwd = 2, col = 'black', 
     ylim = c(0, 1))
lines(0:max_age, downcurve[i, ], lwd = 2, col = 'green')
# some point
abline(h = 0.85, col = 'blue')
abline(v = 2, col = 'blue')
# another point
abline(h = 0.5, col = 'purple')
abline(v = 1, col = 'purple')
# abline(v = 6, col = 'purple')
# where it hits proportion = 1.0
abline(v = 5, col = 'orange')
# abline(h = 0.85, col = 'orange')
# final value
# abline(h = 0.43, col = 'red')