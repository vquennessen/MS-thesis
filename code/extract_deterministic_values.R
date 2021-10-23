## extract deterministic values to overlay on sampling panels

# load any necessary libraries
# library(plyr)
# remotes::install_github('vquennessen/densityratio')
# library(densityratio)

# running on cluster or personal machine?
cluster = TRUE

# determine num_sims based on data folder
num_sims <- 1

# species to compare
species_list <- c('CR_OR_2015', 'BR_OR_2015', 'LING_OW_2017', 'CAB_OR_2019')
Names <- c('Canary Rockfish', 'Black Rockfish', 'Lingcod', 'Cabezon')

# set variables
A = 5
MPA = 3
Time1 = 50
Time2 = 20
Final_DRs1 <- c(0.6, 0.7, 0.8, 0.9)
Final_DRs2 <- c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
Control_rules = c(1:6)
types <- c('Static', 'Transient')
metrics <- c('Biomass', 'Yield', 'Effort')
# , 'Density.Ratio')

nY  <- Time2 + 1
nC  <- length(Control_rules)
nS  <- length(species_list)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFt <- 3*nF1 + nF2
nTy <- length(types)
nM  <- length(metrics)

# set base data frame
DF <- data.frame(Species = c(rep(Names[1], each = nTy*nF1*nY), 
                             rep(Names[2], each = nTy*nF1*nY), 
                             rep(Names[3], each = nTy*nF1*nY), 
                             rep(Names[4], each = nTy*nF2*nY)),
                 Type = c(rep(rep(types, each = nF1*nY), 3), 
                          rep(types, each = nF2*nY)),
                 FDR = c(rep(rep(Final_DRs1, times = nTy, each = nY), 3),
                         rep(Final_DRs2, times = nTy, each = nY)),
                 Year = rep(0:Time2, times = nTy*nFt),
                 Value = rep(NA, nTy*nFt*nY), 
                 Lower = rep(NA, nTy*nFt*nY), 
                 Upper = rep(NA, nTy*nFt*nY))

# copy base dataframes to metric specific ones
BIOMASS <- DF
YIELD   <- DF
EFFORT  <- DF
# DR      <- DF

for (ss in 1:length(species_list)) {
  
  if (cluster == TRUE) {
    # load biomass, yield, and effort files
    load(paste('~/Documents/MS-thesis/data/None/', species_list[ss], 
               '/1_biomass.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/None/', species_list[ss], 
               '/1_yield.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/None/', species_list[ss], 
               '/1_effort.Rda', sep = ''))
    # load(paste('~/Documents/MS-thesis/data/None/', species_list[ss], 
    #            '/1_DR.Rda', sep = ''))
  } else {
    # load biomass, yield, and effort files
    load(paste('~/Projects/MS-thesis/data/None/', species_list[ss], 
               '/1_biomass.Rda', sep = ''))
    load(paste('~/Projects/MS-thesis/data/None/', species_list[ss], 
               '/1_yield.Rda', sep = ''))
    load(paste('~/Projects/MS-thesis/data/None/', species_list[ss], 
               '/1_effort.Rda', sep = ''))
    # load(paste('~/Projects/MS-thesis/data/None/', species_list[ss], 
    #            '/1_DR.Rda', sep = ''))
  }
  
  
  # pull out sample sims as sums across all areas for particular years
  B_sample  <- colSums(sims_biomass[ , , 1:2, , ]) 
  Y_sample  <- sims_yield[ , 1:2, , ]
  E_sample  <- sims_effort[ , 1:2, , ]
  # DR_sample <- sims_DR[ , 1:2, , ]
  
  # set nF based on species
  nF <- ifelse(ss == 4, nF2, nF1)  
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, nY*nTy*nF), c(nY, nTy, nF))
  Rel_yield   <- array(rep(0, nY*nTy*nF), c(nY, nTy, nF))
  Rel_effort  <- array(rep(0, nY*nTy*nF), c(nY, nTy, nF))
  
  # calculate relative arrays after reserve implementation
  for (ty in 1:nTy) {
    for (fdr in 1:nF) {
        Rel_biomass[, ty, fdr] <- B_sample[1:nY, ty, fdr] / B_sample[1, ty, fdr]
        Rel_yield[, ty, fdr] <- Y_sample[1:nY, ty, fdr] / Y_sample[1, ty, fdr]
        Rel_effort[, ty, fdr] <- E_sample[Time1:(Time1 + Time2), ty, fdr] / 
          E_sample[1, ty, fdr]
    }
  }
  
  ##### fill in data frames with median and quantile values #####
  for (ty in 1:nTy) {
    for (fdr in 1:nF) {
      for (y in 1:nY) {
        
        index <- (ss-1)*nTy*nF1*nY + (ty - 1)*nF*nY + (fdr - 1)*nY + y
        
        # relative biomass
        BIOMASS$Value[index] <- Rel_biomass[y, ty, fdr]
        BIOMASS$Lower[index] <- Rel_biomass[y, ty, fdr]
        BIOMASS$Upper[index] <- Rel_biomass[y, ty, fdr]
        
        # relative yield
        YIELD$Value[index] <- Rel_yield[y, ty, fdr]
        YIELD$Lower[index] <- Rel_yield[y, ty, fdr]
        YIELD$Upper[index] <- Rel_yield[y, ty, fdr]
        
        # relative effort
        EFFORT$Value[index] <- Rel_effort[y, ty, fdr]
        EFFORT$Lower[index] <- Rel_effort[y, ty, fdr]
        EFFORT$Upper[index] <- Rel_effort[y, ty, fdr]
        
        # true density ratio
        # DR$Value[index] <- DR_sample[y, ty, fdr]
        # DR$Lower[index] <- DR_sample[y, ty, fdr]
        # DR$Upper[index] <- DR_sample[y, ty, fdr]
        
      }
    }
  }
  
}

##### combine all dataframes into one and write to a .csv file

# add metric columns to dataframes
BIOMASS$Metric <- 'Relative Biomass'
YIELD$Metric   <- 'Relative Yield'
EFFORT$Metric  <- 'Relative Effort'
# DR$Metric      <- 'Density Ratio'

# combine to a single dataframe
# dataframe <- rbind(BIOMASS, YIELD, EFFORT, DR)
DF <- rbind(BIOMASS, YIELD, EFFORT)

# add Source column
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

# rearrange columns to match plot_stochastic dataframe
dataframe <- dataframe[, c('Species', 'Scenario', 'Metric', 'Type', 'FDR', 
                           'Year', 'Value', 'Lower', 'Upper', 'Source')]

# write dataframe to a .csv file
write.csv(x = dataframe, file = 'deterministic_values.csv', row.names = FALSE)
