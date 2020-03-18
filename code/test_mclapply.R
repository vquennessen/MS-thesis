test.mclapply <- function(num_sims, num_cores) {
  
  library(remotes)
  remotes::install_github('vquennessen/densityratio')
  library(densityratio)
  library(parallel)
  library(MASS)
  source('try1.R')
  
  results <- system.time(mclapply(num_sims, try1, mc.cores = num_cores))
  
  filename <- paste('mclapply_', num_sims, 'sims_', num_cores, 'cores.txt')
  write(results, filename, ncolumns = 1)
  
}

