#### selectivity values

library(densityratio)
library(ggplot2)
library(patchwork)
library(egg)

###### canary rockfish #####

Max_age  <- 84                          # maximum age
Rec_age  <- 3                           # age at recruitment
A1       <- 1;         L1   <- 9.05;    # growth parameters (f)
A2       <- 20;        L2   <- 60.05;
K        <- 0.129
Fleets   <- c('trawl', 'non-trawl',     # names of fleets
              'rec', 'hake', 'research')
Alpha    <- c(0.3, 0.6, 1, 1, 1)        # slope of upcurve per fleet
Beta     <- c(1, 0, 1, 1, 0.08)         # slope of downcurve per fleet
F_fin    <- c(0.36, 1, 0.175, 0.65, 0.8)# final selectivity if dome-shaped
A50_up   <- c(5, 5, 4, 8, 1)            # A50 value for upcurve
A50_down <- c(10, 50, 7, 11, 30)        # A50 value for downcurve
Cf       <- c(0.3908, 0.3122, 0.2246,   # fraction of fishery caught / fleet
              0.0295, 0.0429)           #       from upcurve to 1

S1 <- selectivity(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets,A50_up, A50_down, 
                  Alpha, F_fin, Beta, Cf)

ages <- Rec_age:Max_age

DF1 <- data.frame(Age = ages, Selectivity = S1)
DF1$Species <- 'Canary Rockfish'

###### black rockfish #####

Max_age  <- 40                          # maximum age
Rec_age  <- 3                           # age at recruitment
A1       <- 1;        L1 <- 20.32       # growth parameters (f)
A2       <- 40;       L2 <- 49.67
K        <- 0.21
Fleets   <- c('trawl', 'live', 'dead',  # names of fleets
              'ocean', 'shore')
Alpha    <- c(0.325, 0.4, 0.35,
              0.65, 0.425)              # slope of upcurve per fleet
Beta     <- c(0.25, 0.5, 0.4, 1.1, 0.5) # slope of downcurve per fleet
F_fin    <- c(0.325, 0.05, -0.11,
              -0.025, 0.135)            # final selectivity if dome-shaped
A50_up   <- c(7, 5, 5, 5, 3)            # A50 value for upcurve
A50_down <- c(15, 13, 13, 12, 6)        # A50 value for downcurve
Cf       <- c(0.0001, 0.1679, 0.0982,   # fraction of fishery caught / fleet
              0.6979, 0.0358)

S2 <- selectivity(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets,A50_up, A50_down, 
                  Alpha, F_fin, Beta, Cf)

ages <- Rec_age:Max_age

DF2 <- data.frame(Age = ages, Selectivity = S2)
DF2$Species <- 'Black Rockfish'

##### lingcod #####

Max_age  <- 25                          # maximum age
Rec_age  <- 3                           # age at recruitment
A1       <- 1;        L1 <- 17.28;      # growth parameters (f)
A2       <- 20;       L2 <- 120;
K        <- 0.128
Fleets   <- c('trawl', 'fixed_gear',    # names of fleets
              'WArec', 'ORrec')
Alpha    <- c(0.25, 0.25, 0.55, 1)      # slope of upcurve per fleet
Beta     <- c(0.09, 0.3, 0.17, 0.15)    # slope of downcurve per fleet
F_fin    <- c(0.07, 0, 0, 0)            # final select. if dome-shaped
A50_up   <- c(3, 5, 5, 3)               # A50 value for upcurve
A50_down <- c(15, 12, 10, 9)            # A50 value for downcurve
Cf       <- c(0.2872, 0.1379, 0.3253,   # fraction of fishery
              0.2496)

S3 <- selectivity(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets,A50_up, A50_down, 
                  Alpha, F_fin, Beta, Cf)

ages <- Rec_age:Max_age

DF3 <- data.frame(Age = ages, Selectivity = S3)
DF3$Species <- 'Lingcod'

##### cabezon #####

Max_age  <- 20                          # maximum age
Rec_age  <- 4                           # age at recruitment
A1       <- 4;        L1  <- 44.30      # growth parameters (f)
A2       <- 20;       L2  <- 63.35
K <- 0.225
Fleets   <- c('live', 'dead', 'ocean',  # names of fleets
              'shore')
Alpha    <- c(0.4, 0.33, 0.35, 0.9)     # slope of upcurve per fleet
Beta     <- c(0.35, 0, 0, 0.2)          # slope of downcurve per fleet
F_fin    <- c(0.7, 1, 1, 0.07)          # final select. if dome-shaped
A50_up   <- c(3, 4, 2, 1)               # A50 value for upcurve
A50_down <- c(17, 1, 1, 3)              # A50 value for downcurve
Cf       <- c(0.6033, 0.0415, 0.3423,   # fraction of fishery
              0.0130)

S4 <- selectivity(Rec_age, Max_age, A1, L1, A2, L2, K, Fleets,A50_up, A50_down, 
                  Alpha, F_fin, Beta, Cf)

ages <- Rec_age:Max_age

DF4 <- data.frame(Age = ages, Selectivity = S4)
DF4$Species <- 'Cabezon'

# combine dataframes
DF <- rbind(DF1, DF2, DF3, DF4)

g1 <- ggplot(data = DF, aes(x = Age, y = Selectivity)) +
  geom_line() +
  facet_wrap(facets = vars(Species), nrow = 2, ncol = 2, scales = 'free') +
  coord_cartesian(ylim = c(0, 1))

# add tag labels to each facet
g2 <- tag_facet(g1, vjust = 1.15)

# save figure
ggsave('Selectivity_plots.png', g2, 
       path = 'C:/Users/Vic/Box/Quennessen_Thesis/MS Thesis/publication manuscript/figures', 
       width = 6, height = 5)