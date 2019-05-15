#' Victoria Quennessen
#' Calculating Cf for black rockfish fleets

setwd("C:/Users/Vic/Documents/Projects/MS-thesis/code/R/Selectivity parameters")
catch <- read.csv("./black_rockfish_catch.csv", header = T)

year <- catch$Year
n <- length(year)

### OR and CA at once

sport <- catch$Sport_OR + catch$Sport_CA
hook <- catch$Hook_OR + catch$Hook_CA
trawl <- catch$Trawl_OR + catch$Trawl_CA
total <- catch$Total

### Cumulative

c_props <- array(rep(0, 3*n), c(n, 3))

for (i in 1:n) {
  c_props[i, 1] <- sum(sport[i:n])/(sum(total[i:n]))
  c_props[i, 2] <- sum(hook[i:n])/(sum(total[i:n]))
  c_props[i, 3] <- sum(trawl[i:n])/(sum(total[i:n]))
}

c_props <- cbind(year, props)
c_props

### One year at a time

props <- array(rep(0, 3*n), c(n, 3))

for (i in 1:n) {
  props[i, 1] <- sport[i]/total[i]
  props[i, 2] <- hook[i]/total[i]
  props[i, 3] <- trawl[i]/total[i]
}

props <- cbind(year, props)
props

### OR by itself

sport <- catch$Sport_OR 
hook <- catch$Hook_OR 
trawl <- catch$Trawl_OR
total <- sport + hook + trawl

### Cumulative, OR alone

c_props <- array(rep(0, 3*n), c(n, 3))

for (i in 1:n) {
  c_props[i, 1] <- sum(sport[i:n])/(sum(total[i:n]))
  c_props[i, 2] <- sum(hook[i:n])/(sum(total[i:n]))
  c_props[i, 3] <- sum(trawl[i:n])/(sum(total[i:n]))
}

c_props <- cbind(year, c_props)
c_props

### One year at a time, OR alone

props <- array(rep(0, 3*n), c(n, 3))

for (i in 1:n) {
  props[i, 1] <- sport[i]/total[i]
  props[i, 2] <- hook[i]/total[i]
  props[i, 3] <- trawl[i]/total[i]
}

props <- cbind(year, props)
props

### CA by itself

sport <- catch$Sport_CA 
hook <- catch$Hook_CA
trawl <- catch$Trawl_CA
total <- sport + hook + trawl

### Cumulative, CA alone

c_props <- array(rep(0, 3*n), c(n, 3))

for (i in 1:n) {
  c_props[i, 1] <- sum(sport[i:n])/(sum(total[i:n]))
  c_props[i, 2] <- sum(hook[i:n])/(sum(total[i:n]))
  c_props[i, 3] <- sum(trawl[i:n])/(sum(total[i:n]))
}

c_props <- cbind(year, c_props)
c_props

### One year at a time, CA alone

props <- array(rep(0, 3*n), c(n, 3))

for (i in 1:n) {
  props[i, 1] <- sport[i]/total[i]
  props[i, 2] <- hook[i]/total[i]
  props[i, 3] <- trawl[i]/total[i]
}

props <- cbind(year, props)
props
