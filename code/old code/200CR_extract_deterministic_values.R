## extract deterministic values to overlay on sampling panels

# load any necessary libraries
# library(plyr)
remotes::install_github('vquennessen/densityratio')
library(densityratio)

# running on cluster or personal machine?
cluster = TRUE

# determine num_sims based on data folder
num_sims <- 3

# species to compare
species_list <- c('CR_OR_2015')
Names <- c('Canary Rockfish')

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 200
Final_DRs <- c(0.6, 0.7, 0.8, 0.9)
types <- c('Static', 'Transient')
metrics <- c('Relative Biomass', 'Relative Yield', 'Relative Effort')
# , 'Density.Ratio')

nY  <- Time2 + 1
nS  <- length(species_list)
nF  <- length(Final_DRs)
nTy <- length(types)
nM  <- length(metrics)

# set base data frame
DF <- data.frame(Species = rep(Names[1], each = nTy*nF*nY),
                 Type = rep(types, each = nF*nY),
                 FDR = rep(Final_DRs, times = nTy, each = nY),
                 Year = rep(0:Time2, times = nTy*nF),
                 Value = rep(NA, nTy*nF*nY))

# copy base dataframes to metric specific ones
BIOMASS <- DF
YIELD   <- DF
EFFORT  <- DF
DR      <- DF

# load biomass, yield, and effort files
folder <- ifelse(cluster == TRUE, 'Documents', 'Projects')

load(paste('~/', folder, '/MS-thesis/data/None/CR_OR_2015/200 year runs/', 
           '3_biomass.Rda', sep = ''))
load(paste('~/', folder, '/MS-thesis/data/None/CR_OR_2015/200 year runs/', 
           '3_yield.Rda', sep = ''))
load(paste('~/', folder, '/MS-thesis/data/None/CR_OR_2015/200 year runs/', 
           '3_effort.Rda', sep = ''))
load(paste('~/', folder, '/MS-thesis/data/None/CR_OR_2015/200 year runs/', 
           '3_DR.Rda', sep = ''))

# pull out sample sims as sums across all areas for particular years
B_sample  <- colSums(sims_biomass[ , , 1:2, , ]) 
Y_sample  <- sims_yield[ , 1:2, , ]
E_sample  <- sims_effort[ , 1:2, , ]
DR_sample <- sims_DR[ , 1:2, , ]

# initialize relative arrays
Rel_biomass <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
Rel_yield   <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))
Rel_effort  <- array(rep(0, nY*nTy*nF*num_sims), c(nY, nTy, nF, num_sims))

# calculate relative arrays after reserve implementation
for (ty in 1:nTy) {
  for (fdr in 1:nF) {
    for (sim in 1:num_sims) {
      Rel_biomass[, ty, fdr, sim] <- B_sample[, ty, fdr, sim] / 
        B_sample[1, ty, fdr, sim]
      Rel_yield[, ty, fdr, sim] <- Y_sample[, ty, fdr, sim] / 
        Y_sample[1, ty, fdr, sim]
      Rel_effort[, ty, fdr, sim] <- E_sample[Time1:(Time1 + Time2), 
                                             ty, fdr, sim] / 
        E_sample[1, ty, fdr, sim]
    }
  }
}

##### fill in data frames with median and quantile values #####
for (ty in 1:nTy) {
  for (fdr in 1:nF) {
    for (y in 1:nY) {
      
      index <- (ty - 1)*nF*nY + (fdr - 1)*nY + y
      
      # relative biomass
      BIOMASS$Value[index] <- median(Rel_biomass[y, ty, fdr, ])
      
      # relative yield
      YIELD$Value[index] <- median(Rel_yield[y, ty, fdr, ])
      
      # relative effort
      EFFORT$Value[index] <- median(Rel_effort[y, ty, fdr, ])
      
      # true density ratio
      DR$Value[index] <- median(DR_sample[y, ty, fdr, ])
      
    }
  }
}

##### combine all dataframes into one and write to a .csv file

# add metric columns to dataframes
BIOMASS$Metric <- 'Relative Biomass'
YIELD$Metric   <- 'Relative Yield'
EFFORT$Metric  <- 'Relative Effort'
# DR$Metric <- 'Density Ratio'

# combine to a single dataframe
# dataframe <- rbind(BIOMASS, YIELD, EFFORT, DR)
DF <- rbind(BIOMASS, YIELD, EFFORT)

# add Source column
DF$Lower <- NA
DF$Upper <- NA
DF$Source <- 'Deterministic'

# add in scenarios
DF_rec <- DF
DF_rec$Scenario <- 'Recruitment'
DF_samp <- DF
DF_samp$Scenario <- 'Sampling'
DF_both <- DF
DF_both$Scenario <- 'Both'

# merge dataframes again
dataframe <- rbind(DF_rec, DF_samp, DF_both)

# rearrange columns to match plot_sampling_panels dataframe
dataframe <- dataframe[, c('Species', 'Scenario', 'Metric', 'Type', 'FDR', 
                           'Year', 'Value', 'Lower', 'Upper', 'Source')]

# write dataframe to a .csv file
write.csv(x = dataframe, file = '200CR_deterministic_values.csv', 
          row.names = FALSE)
