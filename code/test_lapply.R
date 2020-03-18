test.lapply <- function(num_sims) {
  
  library(remotes)
  remotes::install_github('vquennessen/densityratio')
  library(densityratio)
  library(parallel)
  library(MASS)
  source('try1.R')
  
  results <- system.time(lapply(num_sims, try1))
  
  filename <- paste('lapply_', num_sims, 'sims.txt')
  write(results, filename, ncolumns = 1)
  
}

