initial_size <- function(SAD) {
  
  # find the smallest value in the SAD (typically the plus group)
  Q <- min(SAD)
  
  # find the power you would need to multiply that value by to get a number > 1
  P <- -1*floor(log10(Q))
  
  # set that value as the initial size of the whole population for time = 1, 2
  Init_size <- 1*10^(P/2)
  
  return(Init_size)
  
}