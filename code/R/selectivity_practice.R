# playing around with selectivity

max_age <- 40

start <- 5
peak <- 15
full <- 25
years <- start:max_age

T1 <- 1         # scaling factor
Ffin <- 0.6
b1 <- 0.45      # steepness of ascending curve
b2 <- 7        # inflection point of ascending curve
b3 <- 0.5       # steepness of descending curve
b4 <- 40        # inflection point of descending curve

beta <- array(rep(0, length(years)), c(1, length(years)))

for (age in start:peak){
  beta[age] <- (T1) /( (1 + exp(-1*b1*(age - b2))) )
}

for (age in peak:full){
  beta[age] <- (T1 - Ffin) /( (1 + exp(b3*(age - b4))) )
}

plot(years, beta, ylim = c(0, 1))

max(beta)