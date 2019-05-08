#' Victoria Quennessen
#' Calculating Cf for black rockfish fleets

setwd("C:/Users/Vic/Documents/Projects/MS-thesis/code/R/Selectivity parameters")
catch <- read.csv("./black_rockfish_catch.csv", header = T)

year <- catch$Year
n <- length(year)
sport <- catch$Sport_OR + catch$Sport_CA
hook <- catch$Hook_OR + catch$Hook_CA
trawl <- catch$Trawl_OR + catch$Trawl_CA
total <- catch$Total

c_props <- array(rep(0, 3*n), c(n, 3))

for (i in 1:n) {
  c_props[i, 1] <- sum(sport[i:n])/(sum(total[i:n]))
  c_props[i, 2] <- sum(hook[i:n])/(sum(total[i:n]))
  c_props[i, 3] <- sum(trawl[i:n])/(sum(total[i:n]))
}

c_props <- cbind(year, props)
c_props

props <- array(rep(0, 3*n), c(n, 3))

for (i in 1:n) {
  props[i, 1] <- sport[i]/total[i]
  props[i, 2] <- hook[i]/total[i]
  props[i, 3] <- trawl[i]/total[i]
}

props <- cbind(year, props)
props