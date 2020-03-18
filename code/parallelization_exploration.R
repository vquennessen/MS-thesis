### time for 4 sims normal way: 6962.21 seconds

###### using lapply() to run sims in parallel? #####

rm(list = ls())
setwd('~/Projects/MS-thesis/code')
devtools::install_git('https://github.com/vquennessen/densityratio.git')
library(densityratio)
source('try1.R')

num_sims = 20

system.time(
  lapply(num_sims, try1)
)
 
### time on my laptop (4 sims): 6440.07, 6760.33 seconds
### time on my laptop (20 sims): 18510 seconds
### time on the cluster (4 sims): 3113.632, 2548.881 seconds
### time on the cluster (20 sims): 11522.746,  seconds

###### using mclapply() to run sims in parallel? #####

rm(list = ls())
setwd('~/Projects/MS-thesis/code')
devtools::install_git('https://github.com/vquennessen/densityratio.git')
library(densityratio)
library(parallel)
library(MASS)
source('try1.R')

num_sims = 4

system.time(
  mclapply(num_sims, try1, mc.cores = 4)
)

### cannot be run on Windows
### time on cluster with 4 cores (4 sims): 3487.170, 3362.931
### time on cluster with 12 cores (4 sims): 2742.638, 3176.494

### time on cluster with 4 cores (20 sims): 
### time on cluster with 12 cores (20 sims): 
