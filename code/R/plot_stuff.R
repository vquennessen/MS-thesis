plot_stuff <- function(filepath1, filepath2, filepath3, A, time2, CR, num_sims, 
                       sample_size, PD) {
  
  # setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")  
  
  library(viridis)
  
  load(filepath1)
  load(filepath2)
  load(filepath3)
  
  indices <- sample(1:num_sims, sample_size)
  Y_sample_runs <- sims_yield[, , , indices]
  B_sample_runs <- sims_biomass[, , , indices]
  SSB_sample_runs <- sims_SSB[, , , indices]
  
  # replace NA's in sims_yield with 0s
  Y_sample_runs[is.na(Y_sample_runs)] <- 0
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- Y_lower <- Y_upper <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  
  B_medians <- B_lower <- B_upper <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  
  SSB_medians <- SSB_lower <- SSB_upper <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  
  # extract data from files and plot medians + interquartile ranges
  for (a in 1:A) {
    for (t in 1:time2) {
      for (cr in 1:CR) {

        Y_medians[a, t, cr] <- median(Y_sample_runs[a, t, cr, ])
        B_medians[a, t, cr] <- median(B_sample_runs[a, t, cr, ])
        SSB_medians[a, t, cr] <- median(SSB_sample_runs[a, t, cr, ])
        
        Y_lower[a, t, cr] <- quantile(Y_sample_runs[a, t, cr, ], 0.5 - PD)
        B_lower[a, t, cr] <- quantile(B_sample_runs[a, t, cr, ], 0.5 - PD)
        SSB_lower[a, t, cr] <- quantile(SSB_sample_runs[a, t, cr, ], 0.5 - PD)
         
        Y_upper[a, t, cr] <- quantile(Y_sample_runs[a, t, cr, ], 0.5 + PD)
        B_upper[a, t, cr] <- quantile(B_sample_runs[a, t, cr, ], 0.5 + PD)
        SSB_upper[a, t, cr] <- quantile(SSB_sample_runs[a, t, cr, ], 0.5 + PD)
        
      }
    }
  } 
  
  ###### Plot relative yield + IQR over time after reserve implementation ######  
  
  # use colorblind color palette, viridis
  color <- viridis(CR)
  
  # set plot margins to leave room for legend
  par(mar = c(5.1, 4.1, 4.1, 10.5), xpd = T)
  
  # y-axis limits
  y1 <- 0.2
  y2 <- 1.6
  y_by <- (y2 - y1)/2
  
  # x-axis limits
  x1 <- 0
  x2 <- time2
  x_by <- x2/2
  
  for (a in 1:2) {
    
    area <- ifelse(a == 1, 'far from', 'near')
    title <- sprintf("Relative yield range: %s reserve", area)
    
    # plot the relative biomass
    plot(1, type = 'l',                          # make an empty line graph
         main = title,                           # title of plot
         ylab = 'Yield (metric tons)',              # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                       # set x-axis limits
         ylim = c(y1, y2)
    )
    
    # add a gray dotted line at y = 1
    lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)
    
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
        lines(Y_medians[a, , cr],
              col = color[cr],                  # use pre-defined color palette
              lwd = 2,                     # set line width
              lty = (cr %% 3) + 1)                  # set line type

      polygon(x = c(1:time2, rev(1:time2)), 
              y = c(Y_lower[a, , cr], rev(Y_upper[a, , cr])), 
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
    # add a legend
    legend(x = c(21.25, 28.75), y = c(0.95, 1.65),   # position
           col = color,                           # apply viridis color palette
           lwd = 2,                # apply line thicknesses
           lty = c(2, 3, 1, 2, 3, 1),             # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("B&M Low M", "B&M Correct M", "B&M High M", 
             "Transient Low M", "Transient Correct M", "Transient High M"),           seg.len = 3.5,                        # adjust length of lines
           cex = 0.75)                            # text size
  }
  
  ##### Plot relative biomass + range over time after reserve implementation #####

  # y-axis limits
  yy1 <- 0.5
  yy2 <- 3.8
  yy_by <- (yy2 - yy1)/2

  # x-axis limits
  x1 <- 0
  x2 <- time2
  x_by <- x2/2

  for (a in 1:3) {

    area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
    title <- sprintf("Relative biomass range: %s reserve", area)

    # plot the relative yield
    plot(1, type = 'l',                          # make an empty line graph
         main = title,                           # title of plot
         ylab = 'Relative Biomass (metric tons)',                # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                     # set x-axis limits
         ylim = c(yy1, yy2))

    # add a gray dotted line at y = 1
    lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)
    
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

    for (cr in 1:CR) {
        lines(B_medians[a, , cr],
              col = color[cr],                  # use pre-defined color palette
              lwd = 2,                     # set line width
              lty = c(cr %% 3) + 1)                  # set line type
        
        polygon(x = c(1:time2, rev(1:time2)), 
                y = c(B_lower[a, , cr], rev(B_upper[a, , cr])), 
                col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
    # add a legend
    legend(x = c(21.25, 28.75), y = c(2.35, 3.93),  # position
           col = color,                           # apply viridis color palette
           lwd = 2,                # apply line thicknesses
           lty = c(2, 3, 1, 2, 3, 1),             # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("B&M Low M", "B&M Correct M", "B&M High M", 
             "Transient Low M", "Transient Correct M", "Transient High M"),           seg.len = 3.5,                        # adjust length of lines
           cex = 0.75)                            # text size
  }
  
  ##### Plot relative SSB + range over time after reserve implementation #####

  # y-axis limits
  yyy1 <- 0.5
  yyy2 <- 6.5
  yyy_by <- (yyy2 - yyy1)/2

  # x-axis limits
  x1 <- 0
  x2 <- time2
  x_by <- x2/2

  for (a in 1:3) {

    area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
    title <- sprintf("Relative SSB range: %s reserve", area)

    # plot the relative yield
    plot(1, type = 'l',                          # make an empty line graph
         main = title,                           # title of plot
         ylab = 'Relative SSB (metric tons)',                # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                     # set x-axis limits
         ylim = c(yyy1, yyy2))

    # add a gray dotted line at y = 1
    lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)

    # set specific y-axis
    yyytick <- seq(yyy1, yyy2, by = yyy_by)          # set yaxis tick marks
    axis(side = 2,                               # specify y axis
         at = yyytick,                            # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1)                                # set text horizontal

    # set specific x-axis
    xtick <- seq(x1, x2, by = x_by)          # set x axis tick marks
    axis(side = 1,                               # specify x axis
         at = xtick,                            # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1)                                # set text horizontal

    for (cr in 1:CR) {
        lines(SSB_medians[a, , cr],
              col = color[cr],                  # use pre-defined color palette
              lwd = 2,                     # set line width
              lty = (cr %% 3) + 1)                  # set line type
        
        polygon(x = c(1:time2, rev(1:time2)), 
                y = c(SSB_lower[a, , cr], rev(SSB_upper[a, , cr])), 
                col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
    # add a legend
    legend(x = c(21.25, 28.75), y = c(3.85, 6.7),   # position
           col = color,                           # apply viridis color palette
           lwd = 2,                # apply line thicknesses
           lty = c(2, 3, 1, 2, 3, 1),             # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("B&M Low M", "B&M Correct M", "B&M High M", 
             "Transient Low M", "Transient Correct M", "Transient High M"),           seg.len = 3.5,                        # adjust length of lines
           cex = 0.75)                            # text size
  }

}