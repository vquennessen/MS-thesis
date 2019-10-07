#' Runs base model, based on Babcock & MacCall (2011)

base_model <- function(species, A, time1, time2, allocation, R0, stochasticity, 
                       surveys, fishery_management, fishing, adult_movement) {
  
  # Clear environment and variables
  rm(list = ls())
  
  # load necessary librarys
  library(viridis)
  
  ##### Source functions #######################################################
  
  source("./code/R/parameters.R")
  source("./code/R/length_at_age.R")
  source("./code/R/weight_at_age.R")
  source("./code/R/fraction_mature_at_age.R")
  source("./code/R/selectivity_at_age.R")
  source("./code/R/fishing_mortality.R")
  source("./code/R/epsilon.R")
  source("./code/R/spawning_stock_biomass.R")
  source("./code/R/recruitment.R")
  source("./code/R/pop_dynamics.R")
  source("./code/R/initialize_arrays.R")
  source("./code/R/sampling.R")
  source("./code/R/density_ratio.R")
  source("./code/R/management.R")
  source("./code/R/control_rule.R")
  source("./code/R/Leslie_SAD.R")
  source("./code/R/catch_at_age.R")
  source("./code/R/effort_allocation.R")
  source("./code/R/initial_size.R")
  source("./code/R/vulnerability_to_gear.R")
  source("./code/R/equilibrium_SAD.R")
  source("./code/R/movement.R")
  
  # Set model parameters (fixed)
  CR                   <- 8            # number of control rules
  transects            <- 24           # number of transects per PISCO protocol
                                       #     reserve implementation
  
  ##### Load life history characteristics for species ##########################
  
  par <- parameters(species)
  
  max_age                <- par[[1]]        # maximum age
  M                      <- par[[2]]        # natural mortality
  rec_age                <- par[[3]]        # age at recruitment
  af  <- par[[4]];   bf  <- par[[5]]        # weight at length parameters (f)
  am  <- par[[6]];   bm  <- par[[7]]        # weight at length parameters (m)
  a1f <- par[[8]];  L1f  <- par[[9]]        # growth parameters (f)
  a2f <- par[[10]]; L2f  <- par[[11]] 
  Kf  <- par[[12]]  
  a1m <- par[[13]]; L1m  <- par[[14]]       # growth parameters (m)
  a2m <- par[[15]]; L2m  <- par[[16]]  
  Km                     <- par[[17]]  
  L50                    <- par[[18]]       # length at 50% maturity
  k_mat                  <- par[[19]]       # slope of maturity curve
  ldp                    <- par[[20]]       # larval drift proportion
  h                      <- par[[21]]       # steepness
  phi                    <- par[[22]]       # unfished recruits per spawner
  sigma_R                <- par[[23]]       # recruitment standard deviation
  rho_R                  <- par[[24]]       # recruitment autocorrelation
  AMP                    <- par[[25]]       # adult movement proportion
  D                      <- par[[26]]       # depletion
  Fb                     <- par[[27]]       # fishing mortality to cause D
  r                      <- par[[28]]       # proportion of positive transects 
                                            #       in PISCO monitoring data
  x                      <- par[[29]]       # mean of positive transects
  sp                     <- par[[30]]       # std of positive transects
  c                      <- par[[31]]       # eggs produced per g, intercept
  b                      <- par[[32]]       # eggs produced per g, slope
  
  ####### selectivity parameters #######
  fleets                 <- par[[33]]       # fishery fleet names
  alpha                  <- par[[34]]       # slope for upcurve
  beta                   <- par[[35]]       # slope for downcurve
  start                  <- par[[36]]       # length at initial vulnerability
  F_fin                  <- par[[37]]       # F_fin for fishery, 0 if asymptotic
  L50_up                 <- par[[38]]       # L50 for upcurve
  L50_down               <- par[[39]]       # L50 for downcurve
  cf                     <- par[[40]]       # fraction of fishery caught / fleet
  switch                 <- par[[41]]       # length where selectivity switches 
                                            #       from upcurve to 1
  full                   <- par[[42]]       # length at which downcurve starts
  catch_form             <- par[[43]]       # discrete or continuous catch
  season                 <- par[[44]]       # if catch_formulation = discrete, 
                                            #       time at which fishing occurs:
                                            #       0 at start, 1 at end of year

  ##### Population Dynamics - Non-Time Varying #################################
  
  # Initialize arrays for time-varying dynamics
  IA <- initialize_arrays(A, time1, time2, R0, rec_age, max_age, L1f, L2f, Kf, 
                          a1f, a2f, af, bf, k_mat, Fb, L50, sigma_R, rho_R, 
                          fleets, alpha, beta, start, F_fin, L50_up, L50_down, 
                          cf, switch, full, x, sp, M, CR, phi, catch_form, 
                          season, stochasticity, r, D, transects, h)
  
  timeT            <- IA[[1]]     # total amount of timesteps (years)
  E                <- IA[[2]]     # nominal fishing effort in each area 
  n                <- IA[[3]]     # number of age classes
  L                <- IA[[4]]     # Length at age, dim = 1*age
  W                <- IA[[5]]     # Weight at age, dim = 1*age
  Mat              <- IA[[6]]     # Fraction mature at age, dim = 1*age
  m                <- IA[[7]]     # Age at which fraction mature > 0.5
  S                <- IA[[8]]     # Selectivity at age
  FM               <- IA[[9]]     # Fishing mortality rate, dim = age*area*time
  N                <- IA[[10]]    # Population size, dim = age*area*time
  SSB              <- IA[[11]]    # Spawning stock biomass, dim = area*time
  abundance_all    <- IA[[12]]    # Abundance, dim = area*time
  abundance_mature <- IA[[13]]    # Abundance, dim = area*time
  biomass          <- IA[[14]]    # Biomass, dim = area*time
  Count            <- IA[[15]]    # Species count when sampling, dim = area*time
  nuS              <- IA[[16]]    # Sampling normal variable, dim = area*time*CR
  Eps              <- IA[[17]]    # Epsilon vector, dim = area*time*CR
  catch            <- IA[[18]]    # Catch at age
  yield            <- IA[[19]]    # Yield per area 
  B0               <- IA[[20]]    # Unfished spawning stock biomass
  Delta            <- IA[[21]]    # Constant of proportionality
  Gamma            <- IA[[22]]    # Gamma
  
  ##### Population Dynamics - Time Varying #####################################
  
  for (cr in 1:CR) {
    
    for (t in 3:time1) {
      
      # effort allocation
      E <- effort_allocation(t, cr, allocation, A, E, yield, time1)
      
      # If there is adult movement, add movement
      if (adult_movement == T) { N <- movement(t, cr, N, A, AMP) }
      
      for (a in 1:A) {
        
        # biology
        PD <- pop_dynamics(a, t, cr, rec_age, max_age, n, SSB, N, W, Mat, A, R0, 
                           h, B0, Eps, sigma_R, Fb, E, S, M, FM, m, 
                           abundance_all, abundance_mature, biomass)
        
        SSB                <- PD[[1]]
        FM                 <- PD[[2]]
        N                  <- PD[[3]]
        abundance_all      <- PD[[4]]
        abundance_mature   <- PD[[5]]
        biomass            <- PD[[6]]
        
        # sampling
        if (surveys == T) {
          if (t > (time1 - 3)) {
            Count[a, t, , , cr] <- sampling(a, t, cr, Delta, Gamma, 
                                            abundance_all, abundance_mature, 
                                            transects, x, Count, nuS)
          }
        }
        
        # fishing
        if (fishing == T) {
          catch[, a, t, cr] <- catch_at_age(a, t, cr, FM, M, N, A, Fb, E, catch,
                                            catch_form, season)
          N[, a, t, cr] <- N[, a, t, cr] - catch[, a, t, cr]
          yield[a, t, cr] <- sum(catch[, a, t, cr]*W)
        }
        
      }
      
    }
    
  }
  
  ##### Implement Reserve, and apply control rules #############################
  
  for (cr in 1:CR) {
    
    for (t in (time1 + 1):timeT) {
      
      # effort allocation
      E <- effort_allocation(t, cr, allocation, A, E, yield, time1)
      
      # If there is adult movement, add movement
      if (adult_movement == T) { N <- movement(t, cr, N, A, AMP) }
      
      for (a in 1:A) {
        
        # biology
        PD <- pop_dynamics(a, t, cr, rec_age, max_age, n, SSB, N, W, Mat,
                           A, R0, h, B0, Eps, sigma_R, Fb, E, S, M, FM, m, 
                           abundance_all, abundance_mature, biomass)
        
        SSB                <- PD[[1]]
        FM                 <- PD[[2]]
        N                  <- PD[[3]]
        abundance_all      <- PD[[4]]
        abundance_mature   <- PD[[5]]
        biomass            <- PD[[6]]
        
        # sampling
        if (surveys == T) {
          Count[a, t, , , cr] <- sampling(a, t, cr, Delta, Gamma, abundance_all, 
                                          abundance_mature, transects, x, Count, 
                                          nuS)
        }
        
        # fishing
        if (fishing == T) {
          catch[, a, t, cr] <- catch_at_age(a, t, cr, FM, M, N, A, Fb, E, catch, 
                                            catch_form, season)
          N[, a, t, cr] <- N[, a, t, cr] - catch[, a, t, cr]
          yield[a, t, cr] <- sum(catch[, a, t, cr]*W)
        }
        
      }
      
      # management
      if (fishery_management == T) {
        E <- control_rule(t, cr, E, Count, time1)
      }
      
    }
    
  }
  
  ##### Plot relative biomass over time after reserve implementation ###########
  
  # initialize relative biomass matrix
  rel_biomass <- array(rep(0, A*time2*CR), c(A, time2, CR))
  
  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      rel_biomass[a, , cr] <- biomass[a, (time1 + 1):timeT, cr]/biomass[a, time1, cr]
    }
  }
  
  # use colorblind color palette, viridis
  color <- viridis(CR)
  
  # set plot margins to leave room for legend
  par(mar = c(5.1, 4.1, 4.1, 8.7), xpd = T)
  
  # y-axis limits
  y1 <- 0
  y2 <- 2
  y_by <- (y2 - y1)/2
  
  # x-axis limits
  x1 <- 0
  x2 <- time2
  x_by <- x2/2
  
  for (a in 1:A) {
    title <- sprintf("Relative Biomass per Control Rule: Area %i", a)
    
    # plot the relative biomass
    plot(1, type = 'l',                          # make an empty line graph
         main = title,                           # title of plot
         ylab = 'Relative Biomass',              # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                       # set x-axis limits
         ylim = c(y1, y2)
    )
    
    # set specific y-axis
    ytick <- seq(y1, y2, by = y_by)              # set y axis tick marks
    axis(side = 2,                               # specify y axis
         at = ytick,                             # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1)                                # set text horizontal
    
    # set specific x-axis
    xtick <- seq(x1, x2, by = x_by)              # set x axis tick marks
    axis(side = 1,                               # specify x axis
         at = xtick,                             # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1)                                # set text horizontal    
    
    for (cr in 1:CR) {
      lines(rel_biomass[a, , cr],
            col = color[cr],                     # use pre-defined color palette
            lwd = cr%%2 + 1,                     # set line width
            lty = ceiling(cr/2)                  # set line type
      )
    }
    
    # add a legend
    legend(x = c(54, 62), y = c(y2 + 0.04, y2 - 0.45),   # position
           col = color,                          # apply viridis color palette
           lwd = rep(c(2, 1), 4),                # apply line thicknesses
           lty = rep(1:5, each = 2),             # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("CR 1", "CR 2", "CR 3", "CR 4", "CR 5", "CR 6", "CR 7", "CR 8"),
           seg.len = 3.5,                        # adjust length of lines
           cex = 0.9)                            # text size
  }
  
  ##### Plot relative yield over time after reserve implementation #############
  
  # initialize relative yield matrix
  rel_yield <- array(rep(0, A*time2*CR), c(A, time2, CR))
  
  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      rel_yield[a, , cr] <- yield[a, (time1 + 1):timeT, cr]/yield[a, time1, cr]
    }
  }
  
  # y-axis limits
  yy1 <- -1
  yy2 <- 1 
  yy_by <- (yy2 - yy1)/2
  
  # x-axis limits
  xx1 <- 0
  xx2 <- time2
  xx_by <- xx2/2
  
  for (a in 1:A) {
    title <- sprintf("Relative Yield per Control Rule: Area %i", a)
    
    # plot the relative yield
    plot(1, type = 'l',                          # make an empty line graph
         main = title,                           # title of plot
         ylab = 'Relative Yield',                # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n', 
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(xx1, xx2),                     # set x-axis limits
         ylim = c(yy1, yy2)
    )
    
    
    # set specific y-axis
    yytick <- seq(yy1, yy2, by = yy_by)          # set yaxis tick marks
    axis(side = 2,                               # specify y axis
         at = yytick,                            # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1)                                # set text horizontal
    
    # set specific x-axis
    xxtick <- seq(xx1, xx2, by = xx_by)          # set x axis tick marks
    axis(side = 1,                               # specify x axis
         at = xxtick,                            # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1)                                # set text horizontal    
    
    for (cr in 1:CR) {                           # add one line per control rule
      lines(rel_yield[a, , cr],
            col = color[cr],                     # use pre-defined color palette
            lwd = cr%%2 + 1,                     # set line width
            lty = ceiling(cr/2)                  # set line type
      )
    }
    
    # add a legend
    legend(x = c(54, 62), y = c(yy2 + 0.05, yy2 - 0.55),   # position
           col = color,                          # apply viridis color palette
           lwd = rep(c(2, 1), 4),                # apply line thicknesses
           lty = rep(1:5, each = 2),             # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("CR 1", "CR 2", "CR 3", "CR 4", "CR 5", "CR 6", "CR 7", "CR 8"),
           seg.len = 3.5,                        # adjust length of lines
           cex = 0.9)                            # text size
  }
  
}