catch_at_age <- function(a, t, cr, N, FM) {
  
  catch[ , a, t, cr] <- N[ , a, t, cr]*FM[ , a, t, cr]
  
}