# calculate variability in density ratio calculations

# load packages
library(ggplot2)

##### read in and process data files #####
# read in species data
Species <- read.csv(paste('C:/Users/Vic/Box/Quennessen_Thesis/data/', 
                          'PISCO_taxon_lookup_table.csv', sep = ''), 
                    header = TRUE)

code <- Species$pisco_classcode[min(which(Species$species == 'mystinus'))][1]
# species code for blue rockfish is SMYS

# read in Blue Rockfish data
Fish <- read.csv('C:/Users/Vic/Box/Quennessen_Thesis/data/PISCO_FISH.csv', 
                 header = TRUE)

### notes on PISCO data ###
# site = location - contains either MPA or Reference
# side = ??? - includes CEN (central), DC (down coast), UC (up coast), IN, OUT, 
# OLD_DC, OLD_UC, N, E, S, W
# zone = closer to shore (inner) vs farther (outer)
# level = depth: BOT = bottom, CAN = canopy, MID = midwater
# count = number of individuals found on the transect

# filter out data for other species
BlueRockfish <- subset(Fish, classcode == 'SMYS')

# read in site data
sites_DR <- read.csv(paste('C:/Users/Vic/Box/Quennessen_Thesis/data/',
                           'PISCO subtidal master site table - PISCO sites.csv', 
                           sep = ''))

##### designate sites and samples ##############################################

# sites <- c('Point Lobos SMR', 'Big Creek SMR', 'Vandenberg SMR',
#            'Harris Point SMR', 'South Point SMR')

sites <- c('Point Lobos SMR', 'Harris Point SMR', 'South Point SMR')

num_samples <- 1000

################################################################################

# keep only long-term data
long_term <- subset(sites_DR, SITE_CATEGORY__long_term__Projec == 'LONG-TERM')

# pull out years where both inside and outside zones were sampled
years <- sort(unique(BlueRockfish$year), decreasing = FALSE)

##### calculate density ratios for each site #####

DF <- data.frame(Site = rep(sites, each = length(years)), 
                 Year = rep(years, times = length(sites)), 
                 DR = rep(NA, times = length(sites)*length(years)),
                 SD = rep(NA, times = length(sites)*length(years)))

for (s in 1:length(sites)) {
  
  # make subsets for whole site, and separate ones for inside and outside
  Site <- subset(long_term, MPAGroup == sites[s])
  Site_MPA <- subset(Site, RESERVE == 'IN')
  Site_Ref <- subset(Site, RESERVE == 'OUT')
  MPAs <- unique(as.character(Site_MPA$SITE))
  REFs <- unique(as.character(Site_Ref$SITE))
  IN <- subset(Fish, site %in% MPAs)
  OUT <- subset(Fish, site %in% REFs)
  
  for (y in 1:length(years)) {
    
    # create subset for each year
    YEAR_IN <- subset(IN, year == years[y])
    YEAR_OUT <- subset(OUT, year == years[y])
    
    # check if there are transects inside and outside the reserve
    if (nrow(YEAR_IN) > 0 & nrow(YEAR_OUT) > 0) {
      
      ### separate by transects
      T_in <- split(seq_along(YEAR_IN$transect), with(rle(YEAR_IN$transect), 
                                                      rep(seq_along(values), 
                                                          lengths)))
      T_out <- split(seq_along(YEAR_OUT$transect), with(rle(YEAR_OUT$transect), 
                                                        rep(seq_along(values), 
                                                            lengths)))
      
      # initialize count vector
      T_counts_in <- rep(NA, length(T_in))      
      T_counts_out <- rep(NA, length(T_out))   
      
      # calculate count per transect, 0 if no Blue Rockfish ID
      for (t in 1:length(T_in)) {
        thtr <- YEAR_IN[T_in[[t]], ]
        T_counts_in[t] <- ifelse('SMYS' %in% thtr$classcode,
                                 sum(thtr[which(thtr$classcode == 'SMYS'), 11]), 0)
      }
      
      # calculate count per transect, 0 if no Blue Rockfish ID
      for (t in 1:length(T_out)) {
        thtr <- YEAR_OUT[T_out[[t]], ]
        T_counts_out[t] <- ifelse('SMYS' %in% thtr$classcode,
                                  sum(thtr[which(thtr$classcode == 'SMYS'), 11]), 0)
      }
      
      # calculate final counts and densities    
      # density calculated as sum(count) / # of transects
      Count_in <- sum(T_counts_in)
      Density_in <- Count_in / length(T_in)      
      Count_out <- sum(T_counts_out)
      Density_out <- Count_out / length(T_out)
      
      # calculate simulated variability
      lambda_in <- mean(T_counts_in)
      lambda_out <- mean(T_counts_out)
      DR_simulated <- rep(NA, num_samples)
      
      for (i in 1:num_samples) {
        
        simulated_out <- rpois(n = length(T_in), lambda = lambda_in)
        simulated_in <- rpois(n = length(T_out), lambda = lambda_out)  
        DR_simulated[i] <- mean(simulated_out) / mean(simulated_in)
        
      }
       
      # record observed density ratio and SD in dataframe
      index <- (s - 1)*length(years) + y
      DF$DR[index] <- Density_out / Density_in
      DF$SD[index] <- sd(DR_simulated)
      
    }
    
  }
  
}

# trim useless years off
DF <- na.omit(DF)

# plot DR over time by site
ggplot(DF, aes(x = Year, y = DR)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_errorbar(aes(ymin = DR - SD, ymax = DR + SD), width = 0.25) +
  scale_x_continuous(breaks = seq(2005, 2020, by = 5)) +
  xlim(2004, 2020) +
  facet_grid(cols = vars(Site)) +
  ggtitle('Observed Density Ratio Over Time')
