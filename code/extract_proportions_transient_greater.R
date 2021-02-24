# extract proportions of sims where transient CR led to higher values

# plot relative biomass for static vs. transient DRs for each area

# load any necessary libraries
remotes::install_github('vquennessen/densityratio')
library(densityratio)

###############################################################################
# CHECK THESE EVERY TIME
scenarios <- c('Sampling', 'Recruitment', 'Both')
###############################################################################

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
types <- c('Static', 'Transient')
metrics <- c('Biomass', 'Yield', 'Effort')

nSc <- length(scenarios)
nM  <- length(metrics)
nF1 <- length(Final_DRs1)
nF2 <- length(Final_DRs2)
nFt <- 3*nF1 + nF2
nY  <- Time2 + 1
nTy <- length(types)

base <- data.frame(Species = c(rep(Names[1], each = nSc*nF1*nY), 
                               rep(Names[2], each = nSc*nF1*nY), 
                               rep(Names[3], each = nSc*nF1*nY), 
                               rep(Names[4], each = nSc*nF2*nY)),
                   Scenario = c(rep(rep(scenarios, each = nF1*nY), 3),
                                rep(scenarios, each = nF2*nY)),
                   FDR = c(rep(rep(Final_DRs1, times = nSc, each = nY), 3),
                           rep(Final_DRs2, times = nSc, each = nY)),
                   Year = rep(0:Time2, times = nSc*nFt),
                   Value = rep(NA, nSc*nFt*nY))

# initialize biomass, yield, and effort dataframes
BIO_prop    <- base
YIELD_prop  <- base
EFFORT_prop <- base

# for each species
for (s in 1:length(species_list)) {
  
  # for each scenario
  for (i in 1:length(scenarios)) {
    
    num_sims <- 5000
    
    # load biomass, yield, and effort files
    load(paste('~/Documents/MS-thesis/data/', scenarios[i], '/', 
               species_list[s], '/', num_sims, '_biomass.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/', scenarios[i], '/', 
               species_list[s], '/', num_sims, '_yield.Rda', sep = ''))
    load(paste('~/Documents/MS-thesis/data/', scenarios[i], '/', 
               species_list[s], '/', num_sims, '_effort.Rda', sep = ''))
    
    # set nF value for species 
    nF <- ifelse(s == 4, nF2, nF1)  
    
    ##### relative biomass and median, upper, and lower limits  #####
    
    # pull out sample sims as sums across all areas for particular years
    B_sample   <- colSums(sims_biomass[, , 1:2, , ]) 
    Y_sample <- sims_yield[ , 1:2, , ]
    E_sample <- sims_effort[ , 1:2, , ]
    
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

    for (fdr in 1:nF) {
      for (y in 1:nY) {
        
        index <- (s - 1)*nSc*nF1*nY + (i - 1)*nF*nY + (fdr - 1)*nY + y
        print(index)
        
        # biology dataframe
        BIO_prop$Value[index] <- 
          mean(Rel_biomass[y, 2, fdr, ] > Rel_biomass[y, 1, fdr, ])

        # yield dataframe
        YIELD_prop$Value[index] <- 
          mean(Rel_yield[y, 2, fdr, ] > Rel_yield[y, 1, fdr, ])

        # effort dataframe
        EFFORT_prop$Value[index] <- 
          mean(Rel_effort[y, 2, fdr, ] > Rel_effort[y, 1, fdr, ])

      }
    }
  }
  
}

# add metric columns to dataframes
BIO_prop$Metric    <- 'Biomass'
YIELD_prop$Metric  <- 'Yield'
EFFORT_prop$Metric <- 'Effort'

# combine all dataframes to a single dataframe
DF <- rbind(BIO_prop, YIELD_prop, EFFORT_prop)

# write difference dataframe to csv file #
write.csv(x = DF, file = 'proportion_transient_greater.csv')
