catch_at_age <- function(a, t, N, FM, allocation) {
  
  catch[ , a, t] <- N[, a, t]*FM[, a, t]
  
}