catch_at_age <- function(a, t, N, W, FM, allocation) {
  
  catch[ , a, t] <- N[, a, t]*W*FM[, a, t]
  
}