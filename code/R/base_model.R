#' Runs base model, based on Babcock & MacCall (2011)

base_model <- function(species, A, time1, time2, CR, allocation, R0, 
                       stochasticity, surveys, transects, fishery_management, 
                       fishing, adult_movement, plotting, error, final_DR) {
  
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
  
  # selectivity parameters
  fleets                 <- par[[26]]       # fishery fleet names
  alpha                  <- par[[27]]       # slope for upcurve
  beta                   <- par[[28]]       # slope for downcurve
  F_fin                  <- par[[29]]       # F_fin for fishery, 0 if asymptotic
  L50_up                 <- par[[30]]       # L50 for upcurve
  L50_down               <- par[[31]]       # L50 for downcurve
  cf                     <- par[[32]]       # fraction of fishery caught / fleet
  
  ##### Population Dynamics - Non-Time Varying #################################
  
  # set areas in and out of marine reserves
  areas <- 1:A
  inside <- MPAs
  outside <- areas[-MPAs] 
  
  # Initialize arrays for time-varying dynamics
  IA <- initialize_arrays(A, time1, time2, R0, rec_age, max_age, L1f, L2f, Kf, 
                          a1f, a2f, af, bf, k_mat, Fb, L50, sigma_R, rho_R, 
                          fleets, alpha, beta, start, F_fin, L50_up, L50_down, 
                          cf, switch, full, x, sp, M, CR, phi, stochasticity, 
                          r, D, transects, h, surveys, fishing, error)
  
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
  rel_biomass      <- IA[[23]]    # Relative biomass after reserve implementation
  rel_yield        <- IA[[24]]    # Relative yield after reserve implementation
  rel_SSB          <- IA[[25]]    # Relative SSB after reserve implementation
  nat_mortality    <- IA[[26]]    # Range of potential natural mortality values
  NM               <- IA[[27]]    # Number of potential natural mortality values
  Density_Ratios   <- IA[[28]]    # Density ratios
  
  ##### Population Dynamics - Time Varying #####################################
  
  for (t in 3:time1) {
    
    for (cr in 1:CR) {
      
      for (nm in 1:NM) {
        
        # If there is adult movement, add movement
        if (adult_movement == T) { N <- movement(t, cr, nm, N, A, AMP) }
        
        for (a in 1:A) {
          
          # biology
          PD <- pop_dynamics(a, t, cr, nm, rec_age, max_age, n, SSB, N, W, Mat, 
                             A, R0, h, B0, Eps, sigma_R, Fb, E, S, NM, FM, m, 
                             abundance_all, abundance_mature, biomass, fishing, 
                             nat_mortality)
          
          FM[, a, t, cr, nm]               <- PD[[1]]
          N[, a, t, cr, nm]                <- PD[[2]]
          abundance_all[a, t, cr, nm]      <- PD[[3]]
          abundance_mature[a, t, cr, nm]   <- PD[[4]]
          biomass[a, t, cr, nm]            <- PD[[5]]
          SSB[a, t, cr, nm]                <- PD[[6]]
          
          # sampling
          if (surveys == T) {
              Count[a, t, , , cr, nm] <- sampling(a, t, cr, nm, Delta, Gamma, 
                                                  abundance_all, abundance_mature, 
                                                  transects, x, Count, nuS)
          }
          
          # fishing
          if (fishing == T) {
            catch[, a, t, cr, nm] <- catch_at_age(a, t, cr, nm, FM, nat_mortality, 
                                                  N, A, Fb, E, catch)
            yield[a, t, cr, nm] <- sum(catch[, a, t, cr, nm]*W)
          }
          
        }
      }
      
      # calculate true density ratio
      outside_density <- sum(abundance_mature[c(1, 2, 4, 5), t, cr, 2]) / 4
      inside_density <- sum(abundance_mature[c(3), t, cr, 2]) / 1
      Density_Ratios[t, cr] <- outside_density / inside_density
      
    }
  }
  
  ##### Implement Reserve, and apply control rules #############################
  
  for (t in (time1 + 1):timeT) {
    
    for (cr in 1:CR) {
      
      for (nm in 1:NM) {
        
        # effort allocation
        E <- effort_allocation(t, cr, nm, allocation, A, E, yield, time1)
        
        # If there is adult movement, add movement
        if (adult_movement == T) { N <- movement(t, cr, nm, N, A, AMP) }
        
        for (a in 1:A) {
          
          # biology
          PD <- pop_dynamics(a, t, cr, nm, rec_age, max_age, n, SSB, N, W, Mat, 
                             A, R0, h, B0, Eps, sigma_R, Fb, E, S, NM, FM, m, 
                             abundance_all, abundance_mature, biomass, fishing, 
                             nat_mortality)
          
          FM[, a, t, cr, nm]               <- PD[[1]]
          N[, a, t, cr, nm]                <- PD[[2]]
          abundance_all[a, t, cr, nm]      <- PD[[3]]
          abundance_mature[a, t, cr, nm]   <- PD[[4]]
          biomass[a, t, cr, nm]            <- PD[[5]]
          SSB[a, t, cr, nm]                <- PD[[6]]
          
          # sampling
          if (surveys == T) {
            Count[a, t, , , cr, nm] <- sampling(a, t, cr, nm, Delta, Gamma, 
                                            abundance_all, abundance_mature, 
                                            transects, x, Count, nuS)
          }
          
          # fishing
          if (fishing == T) {
            catch[, a, t, cr, nm] <- catch_at_age(a, t, cr, nm, FM, nat_mortality, 
                                              N, A, Fb, E, catch)
            yield[a, t, cr, nm] <- sum(catch[, a, t, cr, nm]*W)
          }
          
        }
        
        # management
        if (fishery_management == T) {
          E <- control_rule(t, cr, nm, A, E, Count, time1, time2, transects, 
                            nat_mortality, final_DR, inside, outside)
        }
        
      }
      
      # calculate true density ratio
      outside_density <- sum(abundance_mature[c(1, 2, 4, 5), t, cr, 2]) / 4
      inside_density <- sum(abundance_mature[c(3), t, cr, 2]) / 1
      Density_Ratios[t, cr] <- outside_density / inside_density
      
    }
    
  }
  
  # Troubleshooting plots
  # plot(1:timeT, N[1, 1, 1:timeT, 1], type = 'l', col = 'green', ylim = c(0, 2e4))
  # for (x in 2:(n - 1)) {
  #   lines(1:timeT, N[x, 1, 1:timeT, 1], col = 'red')
  # }
  # lines(1:timeT, N[n, 1, 1:timeT, 1], col = 'blue')
  
  ##### Calculate relative biomass, yield, and SSB ############################
  
  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      rel_biomass[a, , cr] <- biomass[a, time1:timeT, cr, 2]/biomass[a, time1, cr, 2]
    }
  }
  
  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      rel_yield[a, , cr] <- yield[a, time1:timeT, cr, 2]/yield[a, time1, cr, 2]
    }
  }
  
  # calculate relative biomass since reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      rel_SSB[a, , cr] <- SSB[a, time1:timeT, cr, 2]/SSB[a, time1, cr, 2]
    }
  }
  
  ##### Plotting ######
  
  if (plotting == T) {
  
    # use red-blue color palette
    palette <- colorRampPalette(c('red', 'blue'))
    color <- palette(CR)
    
    # set line types - solid for correct M, dashed for high M, dotted for low M
    line_type <- c(2, 1, 3, 2, 1, 3)
    
    # set layout matrix for all plots
    layout_m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = T)
    
    # set legend title and text and position
    legend_title <- expression(bold('Control Rule'))
    legend_text  <- c("\n Static \n Low M", "\n Static \n Correct M", 
                     "\n Static \n High M", "\n Transient \n Low M", 
                     "\n Transient \n Correct M", "\n Transient \n High M")
    position <- 'left'
    
    # transient DR for population with correct M
    x_DR <- time2:timeT
    y_DR <- transient_DR(nat_mortality, x_DR, time1, time2, final_DR, nm = 2)
    
    ##### Plot relative biomass over time after reserve implementation #############
    
    # y-axis limits
    y1 <- 0
    y2 <- 4
    y_by <- (y2 - y1)/2
    
    # x-axis limits
    x1 <- 0
    x2 <- time2
    x_by <- x2/4
    
    y1_dr <- 0
    y2_dr <- 2
    by_dr <- y2_dr/2
    
    for (a in 1:3) {

      # set plotting layout
      layout(mat = layout_m,
             widths = c(2, 0.4),                 # Widths of the 2 columns
             heights = c(4, 2))                  # Heights of the 2 rows

      area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
      title <- sprintf("Relative biomass: %s reserve", area)
      
      # plot the relative biomass
      par(mar = c(0.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = title,                         # title of plot
           ylab = 'Relative Biomass',            # axis labels
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(x1, x2),                     # set x-axis limits
           ylim = c(y1, y2), 
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      ytick <- seq(y1, y2, by = y_by)            # set y axis tick marks
      axis(side = 2,                             # specify y axis
           at = ytick,                           # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      for (cr in 1:CR) {
        lines(rel_biomass[a, , cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)               # set line type
      }

      # add a gray dotted line at y = 1
      lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)

      # plot the density ratios over time
      par(mar = c(4.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = 'Density Ratios Over Time',    # title of plot
           ylab = 'Density Ratio',               # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(0, time2),                   # set x-axis limits
           ylim = c(0, y2_dr), 
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      dr_ytick <- seq(y1_dr, y2_dr, by_dr)       # set y axis tick marks
      axis(side = 2,                             # specify y axis
           at = dr_ytick,                        # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      # set specific x-axis
      xtick <- seq(x1, x2, by = x_by)            # set x axis tick marks
      axis(side = 1,                             # specify x axis
           at = xtick,                           # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal
      
      for (cr in 1:CR) {
        lines(Density_Ratios[time1:timeT, cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)}              # set line type
      
      # add a legend
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F, # position
             col = color,                        # apply color palette
             lwd = 2,                            # apply line thicknesses
             lty = line_type,                    # apply line patterns
             title = legend_title,               # add legend title
             legend = legend_text,               # add legend labels
             seg.len = 3,                        # adjust length of lines
             cex = 1.1,                          # adjust legend text size
             bty = 'n') 
      
    }

    ##### Plot relative yield over time after reserve implementation ###########

    # y-axis limits
    yy1 <- 0
    yy2 <- 3
    yy_by <- (yy2 - yy1)/2

    for (a in 1:2) {

      # set plotting layout
      layout(mat = layout_m,
             widths = c(2, 0.4),                 # Widths of the 2 columns
             heights = c(4, 2))                  # Heights of the 2 rows
      
      area <- ifelse(a == 1, 'far from', 'near')
      title <- sprintf("Relative yield: %s reserve", area)
      
      # plot the relative yield
      par(mar = c(0.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = title,                         # title of plot
           ylab = 'Relative Yield',              # axis labels
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(x1, x2),                     # set x-axis limits
           ylim = c(yy1, yy2), 
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      yytick <- seq(yy1, yy2, by = yy_by)        # set yaxis tick marks
      axis(side = 2,                             # specify y axis
           at = yytick,                          # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      for (cr in 1:CR) {
        lines(rel_yield[a, , cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)               # set line type
      }
      
      # add a gray dotted line at y = 1
      lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)

      # plot the density ratio over time
      par(mar = c(4.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = 'Density Ratios Over Time',    # title of plot
           ylab = 'Density Ratio',               # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(0, time2),                   # set x-axis limits
           ylim = c(0, y2_dr), 
           cex.lab = 1.5, cex.main = 1.5)

      # set specific y-axis
      dr_ytick <- seq(y1_dr, y2_dr, by_dr)       # set y axis tick marks
      axis(side = 2,                             # specify y axis
           at = dr_ytick,                        # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal

      # set specific x-axis
      xtick <- seq(x1, x2, by = x_by)            # set x axis tick marks
      axis(side = 1,                             # specify x axis
           at = xtick,                           # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal
      
      for (cr in 1:CR) {
        lines(Density_Ratios[time1:timeT, cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)}              # set line type
      
      # add a legend
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F, # position
             col = color,                        # apply color palette
             lwd = 2,                            # apply line thicknesses
             lty = line_type,                    # apply line patterns
             title = legend_title,               # add legend title
             legend = legend_text,               # add legend labels
             seg.len = 3,                        # adjust length of lines
             cex = 1.1,                          # adjust legend text size
             bty = 'n') 
      
    }

    ###### Plot relative SSB over time after reserve implementation ##############
    
    # y-axis limits
    yyy1 <- 0
    yyy2 <- 2 
    yyy_by <- (yyy2 - yyy1)/2
    
    for (a in 1:3) {
      
      # set plotting layout
      layout(mat = layout_m,
             widths = c(2, 0.4),                 # Widths of the 2 columns
             heights = c(4, 2))                  # Heights of the 2 rows
      
      area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
      title <- sprintf("Relative biomass: %s reserve", area)
      
      # plot the relative SSB
      par(mar = c(0.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = title,                         # title of plot
           ylab = 'Relative SSB',                # axis labels
           xaxt = 'n', 
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(x1, x2),                     # set x-axis limits
           ylim = c(yyy1, yyy2), 
           cex.lab = 1.5, cex.main = 1.5)
      
      # set specific y-axis
      yyytick <- seq(yyy1, yyy2, by = yyy_by)    # set yaxis tick marks
      axis(side = 2,                             # specify y axis
           at = yyytick,                         # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal
      
      for (cr in 1:CR) {
        lines(rel_SSB[a, , cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)               # set line type
      }
      
      # add a gray dotted line at y = 1
      lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)
      
      # plot the density ratio over time
      par(mar = c(4.1, 4.5, 3.1, 0.1))
      plot(1, type = 'l',                        # make an empty line graph
           main = 'Density Ratios Over Time',    # title of plot
           ylab = 'Density Ratio',               # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                           # get rid of y-axis
           xlim = c(0, time2),                   # set x-axis limits
           ylim = c(0, y2_dr),
           cex.lab = 1.5, cex.main = 1.5
      )
      
      # set specific y-axis
      dr_ytick <- seq(y1_dr, y2_dr, by_dr)       # set y axis tick marks
      axis(side = 2,                             # specify y axis
           at = dr_ytick,                        # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal
      
      # set specific x-axis
      dr_xtick <- seq(0, time2, by = time2/4)    # set x axis tick marks
      axis(side = 1,                             # specify x axis
           at = dr_xtick,                        # apply tick marks
           labels = T,                           # apply appropriate labels
           las = 1)                              # set text horizontal    
      
      for (cr in 1:CR) {
        lines(Density_Ratios[time1:timeT, cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = (cr %% 3) + 1)}              # set line type
      
      
      # add a legend
      par(mar = c(0.1, 0.1, 0.1, 0.1))
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F, # position
             col = color,                        # apply color palette
             lwd = 2,                            # apply line thicknesses
             lty = line_type,                    # apply line patterns
             title = legend_title,               # add legend title
             legend = legend_text,               # add legend labels
             seg.len = 3,                        # adjust length of lines
             cex = 1.1,                          # adjust legend text size
             bty = 'n') 
      
      }
    
  }
  
  output <- list(rel_yield, rel_biomass, rel_SSB, Density_Ratios)
  
  return(output)
  
}