#' Runs base model, based on Babcock & MacCall (2011)

base_model <- function(species, A, time1, time2, CR, allocation, R0, 
                       stochasticity, surveys, transects, fishery_management, 
                       fishing, adult_movement, plotting) {
  
  # load necessary librarys
  if (plotting == T) {library(viridis)}
  
  ##### Source functions #######################################################
  
  source("./parameters.R")
  source("./length_at_age.R")
  source("./weight_at_age.R")
  source("./fraction_mature_at_age.R")
  source("./selectivity_at_age.R")
  source("./fishing_mortality.R")
  source("./epsilon.R")
  source("./spawning_stock_biomass.R")
  source("./recruitment.R")
  source("./pop_dynamics.R")
  source("./initialize_arrays.R")
  source("./sampling.R")
  source("./density_ratio.R")
  source("./management.R")
  source("./control_rule.R")
  source("./Leslie_SAD.R")
  source("./catch_at_age.R")
  source("./effort_allocation.R")
  source("./equilibrium_SAD.R")
  source("./movement.R")
  source("./transient_DR.R")
  
  ##### Load life history characteristics for species ##########################
  
  par <- parameters(species)
  
  max_age                <- par[[1]]        # maximum age
  M                      <- par[[2]]        # natural mortality
  rec_age                <- par[[3]]        # age at recruitment
  af  <- par[[4]];   bf  <- par[[5]]        # weight at length parameters (f)
  a1f <- par[[6]];  L1f  <- par[[7]]        # growth parameters (f)
  a2f <- par[[8]];  L2f  <- par[[9]] 
  Kf  <- par[[10]]  
  L50                    <- par[[11]]       # length at 50% maturity
  k_mat                  <- par[[12]]       # slope of maturity curve
  ldp                    <- par[[13]]       # larval drift proportion
  h                      <- par[[14]]       # steepness
  phi                    <- par[[15]]       # unfished recruits per spawner
  sigma_R                <- par[[16]]       # recruitment standard deviation
  rho_R                  <- par[[17]]       # recruitment autocorrelation
  AMP                    <- par[[18]]       # adult movement proportion
  D                      <- par[[19]]       # depletion
  Fb                     <- par[[20]]       # fishing mortality to cause D
  r                      <- par[[21]]       # proportion of positive transects 
                                            #       in PISCO monitoring data
  x                      <- par[[22]]       # mean of positive transects
  sp                     <- par[[23]]       # std of positive transects
  c                      <- par[[24]]       # eggs produced per g, intercept
  b                      <- par[[25]]       # eggs produced per g, slope
  
  ####### selectivity parameters #######
  fleets                 <- par[[26]]       # fishery fleet names
  alpha                  <- par[[27]]       # slope for upcurve
  beta                   <- par[[28]]       # slope for downcurve
  F_fin                  <- par[[29]]       # F_fin for fishery, 0 if asymptotic
  L50_up                 <- par[[30]]       # L50 for upcurve
  L50_down               <- par[[31]]       # L50 for downcurve
  cf                     <- par[[32]]       # fraction of fishery caught / fleet
  
  ##### Population Dynamics - Non-Time Varying #################################
  
  # Initialize arrays for time-varying dynamics
  IA <- initialize_arrays(A, time1, time2, R0, rec_age, max_age, L1f, L2f, Kf, 
                          a1f, a2f, af, bf, k_mat, Fb, L50, sigma_R, rho_R, 
                          fleets, alpha, beta, start, F_fin, L50_up, L50_down, 
                          cf, switch, full, x, sp, M, CR, phi, stochasticity, r, 
                          D, transects, h)
  
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
      
      # If there is adult movement, add movement
      if (adult_movement == T) { N <- movement(t, cr, N, A, AMP) }
      
      for (a in 1:A) {
        
        # biology
        PD <- pop_dynamics(a, t, cr, rec_age, max_age, n, SSB, N, W, Mat, A, R0, 
                           h, B0, Eps, sigma_R, Fb, E, S, M, FM, m, 
                           abundance_all, abundance_mature, biomass)
        
        SSB                <- PD[[1]]
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
          catch[, a, t, cr] <- catch_at_age(a, t, cr, FM, M, N, A, Fb, E, catch)
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
          catch[, a, t, cr] <- catch_at_age(a, t, cr, FM, M, N, A, Fb, E, catch)
          N[, a, t, cr] <- N[, a, t, cr] - catch[, a, t, cr]
          yield[a, t, cr] <- sum(catch[, a, t, cr]*W)
        }
        
      }
      
      # management
      if (fishery_management == T) {
        E <- control_rule(t, cr, E, Count, time1, time2, transects, M)
      }
      
    }
    
  }
  
  # initialize relative biomass matrix
  rel_biomass <- array(rep(0, A*time2*CR), c(A, time2, CR))
  
  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      rel_biomass[a, , cr] <- biomass[a, (time1 + 1):timeT, cr]/biomass[a, time1, cr]
    }
  }
  
  # initialize relative yield matrix
  rel_yield <- array(rep(0, A*time2*CR), c(A, time2, CR))
  
  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      rel_yield[a, , cr] <- yield[a, (time1 + 1):timeT, cr]/yield[a, time1, cr]
    }
  }
  
  # initialize relative spawning stock biomass matrix
  rel_SSB <- array(rep(0, A*time2*CR), c(A, time2, CR))
  
  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      rel_SSB[a, , cr] <- SSB[a, (time1 + 1):timeT, cr]/SSB[a, time1, cr]
    }
  }
  
  if (plotting == T) {
    
    ##### Plot relative biomass over time after reserve implementation #########
    
    # use colorblind color palette, viridis
    color <- viridis(CR)
    
    # set plot margins to leave room for legend
    par(mar = c(5.1, 4.1, 4.1, 8.7), xpd = T)
    
    # y-axis limits
    y1 <- 0.5
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
              col = color[cr],                  # use pre-defined color palette
              lwd = 2,                     # set line width
              lty = cr)                  # set line type
      }
      
      # add a legend
      legend(x = c(22, 28), y = c(y2 + 0.04, y2 - 0.45),   # position
             col = color,                          # apply viridis color palette
             lwd = 2,                # apply line thicknesses
             lty = 1:CR,             # apply line patterns
             title = 'CR',                         # add legend title and labels
             c("CR 1", "CR 2", "CR 3", "CR 4", "CR 5"),
             seg.len = 3.5,                        # adjust length of lines
             cex = 0.9)                            # text size
    }
    
    ##### Plot relative yield over time after reserve implementation ###########
  
    # y-axis limits
    yy1 <- 0
    yy2 <- 3 
    yy_by <- (yy2 - yy1)/2
    
    for (a in 1:A) {
      title <- sprintf("Relative Yield per Control Rule: Area %i", a)
      
      # plot the relative yield
      plot(1, type = 'l',                          # make an empty line graph
           main = title,                           # title of plot
           ylab = 'Relative Yield',                # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n', 
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(x1, x2),                     # set x-axis limits
           ylim = c(yy1, yy2)
      )
      
      
      # set specific y-axis
      yytick <- seq(yy1, yy2, by = yy_by)          # set yaxis tick marks
      axis(side = 2,                               # specify y axis
           at = yytick,                            # apply tick marks
           labels = T,                             # apply appropriate labels
           las = 1)                                # set text horizontal
      
      # set specific x-axis
      xtick <- seq(x1, x2, by = x_by)          # set x axis tick marks
      axis(side = 1,                               # specify x axis
           at = xtick,                            # apply tick marks
           labels = T,                             # apply appropriate labels
           las = 1)                                # set text horizontal    
      
      for (cr in 1:CR) {                        # add one line per control rule
        lines(rel_yield[a, , cr],
              col = color[cr],                  # use pre-defined color palette
              lwd = 2,                     # set line width
              lty = cr)                  # set line type
      }
      
      # add a legend
      legend(x = c(22, 28), y = c(yy2 + 0.05, yy2 - 0.55),   # position
             col = color,                          # apply viridis color palette
             lwd = 2,                # apply line thicknesses
             lty = 1:CR,             # apply line patterns
             title = 'CR',                         # add legend title and labels
             c("CR 1", "CR 2", "CR 3", "CR 4", "CR 5"),
             seg.len = 3.5,                        # adjust length of lines
             cex = 0.9)                            # text size
    }
    
  }
  
  output <- list(rel_yield, rel_biomass, rel_SSB)
  
  return(output)
  
}