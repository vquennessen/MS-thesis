################################################################################
##### Calculate alpha values
################################################################################

# clear environment
rm(list = ls())

# source required functions
source("./code/R/parameters.R")
source("./code/R/length_at_age.R")

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
source("./code/R/parameters.R")
source("./code/R/length_at_age.R")

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
source("./code/R/parameters.R")
source("./code/R/length_at_age.R")

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
alpha <- c(0.8473, 0.6931, 2.1972, 1.3863, 2.1972)
Ffin <- c(0.6, 0.275, 0.2, 0.015, 0.13)
L50downA <- c(15, 13, 14, 12, 6.5)
L50downL <- L[L50downA - rec_age + 1]
beta <- c(-0.1703, 0.5229, 0.6486, 0.5593, 0.1231)

# initialize upcurves and downcurves
f <- length(fleet)
l <- length(L)
upcurve <- array(rep(NA, f * l), c(f, l))
downcurve <- array(rep(NA, f * l), c(f, l))

for (i in 1:l) {
  upcurve[i, ] <- 1 / (1 + exp(-1 * alpha[i] * (L - L50upL[i])))
  downcurve[i, ] <- 1 - (1 - Ffin[i]) / (1 + exp(-1 * beta[i] * (L - L50downL[i])))
}

plot(age, upcurve[1, ], type = 'l', lwd = 2, col = 'darkblue')
