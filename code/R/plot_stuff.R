# plot_stuff <- function(filepath1, filepath2, filepath3, A, time2, CR, num_sims, 
#                        sample_size, PD) {
  
  setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")  
  
  library(viridis)
  
  load(filepath1)
  load(filepath2)
  load(filepath3)
  
  A <- 5
  CR <- 5
  time2 <- 20
  num_sims <- 1e4
  sample_size <- 5e2
  PD <- 0.2
  
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
  for (i in 1:A) {
    for (j in 1:time2) {
      for (k in 1:CR) {
        
        Y_medians[i, j, k] <- median(Y_sample_runs[i, j, k, ])
        B_medians[i, j, k] <- median(B_sample_runs[i, j, k, ])
        SSB_medians[i, j, k] <- median(SSB_sample_runs[i, j, k, ])
        
        Y_lower[i, j, k] <- quantile(Y_sample_runs[i, j, k, ], 0.5 - PD)
        B_lower[i, j, k] <- quantile(B_sample_runs[i, j, k, ], 0.5 - PD)
        SSB_lower[i, j, k] <- quantile(SSB_sample_runs[i, j, k, ], 0.5 - PD)
         
        Y_upper[i, j, k] <- quantile(Y_sample_runs[i, j, k, ], 0.5 + PD)
        B_upper[i, j, k] <- quantile(B_sample_runs[i, j, k, ], 0.5 + PD)
        SSB_upper[i, j, k] <- quantile(SSB_sample_runs[i, j, k, ], 0.5 + PD)
        
      }
    }
  } 
  
  ###### Plot relative yield + IQR over time after reserve implementation ######  
  
  # use colorblind color palette, viridis
  color <- viridis(CR)
  
  # set plot margins to leave room for legend
  par(mar = c(5.1, 4.1, 4.1, 8.7), xpd = T)
  
  # y-axis limits
  y1 <- 0
  y2 <- 1.4
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
            col = color[cr],                     # use pre-defined color palette
            lwd = 2,                     # set line width
            lty = cr)                  # set line type
      
      polygon(x = c(1:time2, rev(1:time2)), 
              y = c(Y_lower[a, , cr], rev(Y_upper[a, , cr])), 
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
    # add a legend
    legend(x = c(x2 + 1.5, x2 + 11), y = c(y2, y2 - (y2-y1)/2),   # position
           col = color,                          # apply viridis color palette
           lwd = 2,                # apply line thicknesses
           lty = 1:CR,             # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("B&M Best", "B&M Worst", "Low M", "Correct M", "High M"),
           seg.len = 3.5,                        # adjust length of lines
           cex = 0.9)                            # text size
  }
  
  ##### Plot relative biomass + range over time after reserve implementation #####

  # y-axis limits
  yy1 <- 0
  yy2 <- 4
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

    for (cr in 1:CR) {                           # add one line per control rule
      lines(B_medians[a, , cr],
            col = color[cr],                     # use pre-defined color palette
            lwd = 2,                     # set line width
            lty = cr)                 # set line type

      polygon(x = c(1:time2, rev(1:time2)),
              y = c(B_lower[a, , cr], rev(B_upper[a, , cr])),
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }

    # add a legend
    legend(x = c(x2 + 1.5, x2 + 11), y = c(yy2 + 0.04, yy2 - (yy2-yy1)/2.1),   # position
           col = color,                          # apply viridis color palette
           lwd = 2,                              # apply line thicknesses
           lty = 1:CR,                            # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("B&M Best", "B&M Worst", "Low M", "Correct M", "High M"),
           seg.len = 3.5,                        # adjust length of lines
           cex = 0.9)                            # text size
  }
  
  ##### Plot relative SSB + range over time after reserve implementation #####

  # y-axis limits
  yyy1 <- 0
  yyy2 <- 7
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
         ylim = c(yyy1, yyy2)
    )


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

    for (cr in 1:CR) {                           # add one line per control rule
      lines(SSB_medians[a, , cr],
            col = color[cr],                     # use pre-defined color palette
            lwd = 2,                     # set line width
            lty = cr)                 # set line type

      polygon(x = c(1:time2, rev(1:time2)),
              y = c(SSB_lower[a, , cr], rev(SSB_upper[a, , cr])),
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }

    # add a legend
    legend(x = c(x2 + 1.5, x2 + 11), y = c(yyy2 + 0.04, yyy2 - (yyy2-yyy1)/2.1),   # position
           col = color,                          # apply viridis color palette
           lwd = 2,                              # apply line thicknesses
           lty = 1:CR,                            # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("B&M Best", "B&M Worst", "Low M", "Correct M", "High M"),
           seg.len = 3,                        # adjust length of lines
           cex = 0.9)                            # text size
  }

# }
