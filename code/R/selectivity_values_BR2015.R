################################################################################
##### Calculate alpha values
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'black rockfish 2015'

fleet <- c('trawl', 'live', 'dead', 'recO', 'recS')
L50up <- c(7, 5, 5, 5, 3)
Sa <- c(0.3, 0.8, 0.1, 0.2, 0.1)
La <- c(6, 7, 4, 4, 2)

# load species parameters
par <- parameters(species)

max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
rec_age                <- par[[3]]        # age at recruitment
af  <- par[[4]];   bf  <- par[[5]]        # weight at length parameters (f)
am  <- par[[6]];   bm  <- par[[7]]        # weight at length parameters (m)
a1f <- par[[8]];  L1f  <- par[[9]]        # growth parameters (f)
a2f <- par[[10]]; L2f  <- par[[11]] 
Kf  <- par[[12]]  
a1m <- par[[13]]; L1m  <- par[[14]]       # growth parameters (m)
a2m <- par[[15]]; L2m  <- par[[16]]  
Km                     <- par[[17]]  
L50                    <- par[[18]]       # length at 50% maturity
k_mat                  <- par[[19]]       # slope of maturity curve

# Calculated values
age <- rec_age:max_age                          # applicable ages
n <- length(age)                                # number of age bins
L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f) # length at age

alpha <- (-1 * log(1/Sa - 1)) / (La - L50up)

################################################################################
##### Calculate beta values
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'black rockfish 2015'

fleet <- c('trawl', 'live', 'dead', 'recO', 'recS')
Ffin <- c(0.6, 0.275, 0.2, 0.015, 0.13)
L50down <- c(15, 13, 14, 12, 6.5)
Sa <- c(0.75, 0.4, 0.3, 0.17, 0.3)
La <- c(12, 16, 17, 15, 18)

# load species parameters
par <- parameters(species)

max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
rec_age                <- par[[3]]        # age at recruitment
af  <- par[[4]];   bf  <- par[[5]]        # weight at length parameters (f)
am  <- par[[6]];   bm  <- par[[7]]        # weight at length parameters (m)
a1f <- par[[8]];  L1f  <- par[[9]]        # growth parameters (f)
a2f <- par[[10]]; L2f  <- par[[11]] 
Kf  <- par[[12]]  
a1m <- par[[13]]; L1m  <- par[[14]]       # growth parameters (m)
a2m <- par[[15]]; L2m  <- par[[16]]  
Km                     <- par[[17]]  
L50                    <- par[[18]]       # length at 50% maturity
k_mat                  <- par[[19]]       # slope of maturity curve

# Calculated values
age <- rec_age:max_age                          # applicable ages
n <- length(age)                                # number of age bins
L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f) # length at age

beta <- (-1 * log((1 - Ffin)/(1 - Sa) - 1)) / (La - L50down)

################################################################################
##### try to graph selectivity
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./parameters.R")
source("./length_at_age.R")

species <- 'black rockfish 2015'

# load species parameters
par <- parameters(species)

max_age                <- par[[1]]        # maximum age
M                      <- par[[2]]        # natural mortality
rec_age                <- par[[3]]        # age at recruitment
af  <- par[[4]];   bf  <- par[[5]]        # weight at length parameters (f)
am  <- par[[6]];   bm  <- par[[7]]        # weight at length parameters (m)
a1f <- par[[8]];  L1f  <- par[[9]]        # growth parameters (f)
a2f <- par[[10]]; L2f  <- par[[11]] 
Kf  <- par[[12]]  
a1m <- par[[13]]; L1m  <- par[[14]]       # growth parameters (m)
a2m <- par[[15]]; L2m  <- par[[16]]  
Km                     <- par[[17]]  

# Calculated values
age <- rec_age:max_age                          # applicable ages
n <- length(age)                                # number of age bins
L <- length_at_age(age, L1f, L2f, Kf, a1f, a2f) # length at age

# manual parameters
fleet <- c('trawl', 'live', 'dead', 'recO', 'recS')
L50upA <- c(7, 5, 5, 5, 3)
L50upL <- L[L50upA - rec_age + 1]
alpha <- c(0.325, 0.4, 0.35, 0.65, 0.425)
Ffin <- c(0.325, 0.05, -0.11, -0.025, 0.135)
L50downA <- c(15, 13, 13, 12, 6.5)
L50downL <- L[L50downA - rec_age + 1]
beta <- c(0.25, 0.5, 0.4, 1.1, 0.5)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)
upcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))
downcurve <- array(rep(NA, f * (max_age + 1)), c(f, max_age + 1))
for (i in 1:rec_age) {
  upcurve[, i] <- downcurve[, i] <- 0
}

for (i in 1:f) {
  upcurve[i, age + 1] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50upL[i])))
  downcurve[i, age + 1] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50downL[i])))
}

x <- 5
plot(0:max_age, upcurve[x, ], type = 'l', lwd = 2, col = 'red', 
     ylim = c(0, 1))
lines(0:max_age, downcurve[x, ], lwd = 2, col = 'red')
abline(v = 4.5, col = 'red')
abline(h = 0.875, col = 'red')
abline(h = 0.14, col = 'green')
abline(h = 0.5, col = 'blue')
abline(v = 3, col = 'blue')
abline(v = 6.25, col = 'blue')