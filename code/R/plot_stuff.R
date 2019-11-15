plot_stuff <- function(filepath1, filepath2, filepath3, filepath4, 
                       A, time2, CR, num_sims, sample_size, PD, 
                       plot_individual_runs) {
  
  # setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")  
  
  # library(viridis)
  
  load(filepath1)
  load(filepath2)
  load(filepath3)
  load(filepath4)
  
  indices <- sample(1:num_sims, sample_size)
  Y_sample_runs <- sims_yield[, , , indices]
  B_sample_runs <- sims_biomass[, , , indices]
  SSB_sample_runs <- sims_SSB[, , , indices]
  DR_sample_runs <- sims_DR[, , indices]
  
  # replace NA's in sims_yield with 0s
  Y_sample_runs[is.na(Y_sample_runs)] <- 0
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- Y_lower <- Y_upper <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  B_medians <- B_lower <- B_upper <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  SSB_medians <- SSB_lower <- SSB_upper <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  DR_medians <- array(rep(NA, time2*CR), c(time2, CR))
  
  # extract data from files and plot medians + interquartile ranges
  for (t in 1:time2) {
    for (cr in 1:CR) {
      
      DR_medians[t, cr] <- median(DR_sample_runs[t, cr, ])
      
      for (a in 1:A) {
        
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
  
  # use red-blue color palette
  palette <- colorRampPalette(c('red', 'blue'))
  color <- palette(CR)
  
  # set line types - solid for correct M, dashed for high M, dotted for low M
  line_type <- c(2, 1, 3, 2, 1, 3)
  
  # y-axis limits
  y1 <- 0.3
  y2 <- 1.5
  y_by <- (y2 - y1)/4
  
  # x-axis limits
  x1 <- 0
  x2 <- time2
  x_by <- x2/4
  
  # density ratio y limits
  y1_dr <- 0
  y2_dr <- 3
  by_dr <- y2_dr/2
  
  for (a in 1:2) {
    
    area <- ifelse(a == 1, 'far from', 'near')
    title <- sprintf("Relative yield: %s reserve", area)
    
    # set plotting layout
    m <- matrix(c(1, 2, 3), nrow = 3, ncol = 1, byrow = T)
    layout(mat = m,
           widths = c(2),                      # Widths of the 3 columns
           heights = c(3.75, 2, 1))               # Heights of the 2 rows
    
    # plot the relative biomass
    par(mar = c(0.1, 4.5, 3.1, 0.1))
    plot(1, type = 'l',                          # make an empty line graph
         main = title,                           # title of plot
         ylab = 'Relative Yield',              # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                       # set x-axis limits
         ylim = c(y1, y2), 
         cex.main = 1.75, cex.axis = 3, cex.lab = 1.5)
    
    # add a gray dotted line at y = 1
    lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)
    
    # set specific y-axis
    ytick <- seq(y1, y2, by = y_by)              # set y axis tick marks
    axis(side = 2,                               # specify y axis
         at = ytick,                             # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1)                                # set text horizontal
    
    for (cr in 1:CR) {
      lines(Y_medians[a, , cr],
            col = color[cr],                  # use pre-defined color palette
            lwd = 2,                     # set line width
            lty = line_type[cr])                  # set line type
      
      polygon(x = c(1:time2, rev(1:time2)), 
              y = c(Y_lower[a, , cr], rev(Y_upper[a, , cr])), 
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
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
    dr_ytick <- seq(y1_dr, y2_dr, by_dr)              # set y axis tick marks
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
      lines(DR_medians[, cr],
            col = color[cr],                   # use pre-defined color palette
            lwd = 2,                           # set line width
            lty = (cr %% 3) + 1)}              # set line type
    
    # add a gray dotted line at y = final_target_DR
    lines(0:time2, rep(0.6, time2 + 1), col = 'gray', lty = 3)
    
    # add a legend
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(1, type = 'n', axes = F, xlab = '', ylab = '')
    legend(x = 'top', inset = 0, horiz = T,    # position
           col = color,                              # apply color palette
           lwd = 2,                                  # apply line thicknesses
           lty = line_type,                          # apply line patterns
           title = expression(bold('Control Rule')), # add legend title, labels
           c("Static \n Low M", "Static \n Correct M", "Static \n High M",
             "Transient \n Low M", "Transient \n Correct M", "Transient \n High M"),           
           seg.len = 3,                              # adjust length of lines
           cex = 1.1,                                # text size
           bty = 'n') 
    
  }
  
  ##### Plot relative biomass + range over time after reserve implementation #####
  
  # y-axis limits
  yy1 <- 0.6
  yy2 <- 3.6
  yy_by <- (yy2 - yy1)/4
  
  for (a in 1:3) {
    
    # set plotting layout
    m <- matrix(c(1, 2, 3), nrow = 3, ncol = 1, byrow = T)
    layout(mat = m,
           widths = c(2),                      # Widths of the 3 columns
           heights = c(3.75, 2, 1))               # Heights of the 2 rows
    
    area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
    title <- sprintf("Relative biomass: %s reserve", area)
    
    # plot the relative yield
    par(mar = c(0.1, 4.5, 3.1, 0.1))
    plot(1, type = 'l',                          # make an empty line graph
         main = title,                           # title of plot
         ylab = 'Relative Biomass',                # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                     # set x-axis limits
         ylim = c(yy1, yy2), 
         cex.main = 1.75, cex.axis = 3, cex.lab = 1.5)
    
    # add a gray dotted line at y = 1
    lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)
    
    # set specific y-axis
    yytick <- seq(yy1, yy2, by = yy_by)          # set yaxis tick marks
    axis(side = 2,                               # specify y axis
         at = yytick,                            # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1)                                # set text horizontal
    
    for (cr in 1:CR) {
      lines(B_medians[a, , cr],
            col = color[cr],                  # use pre-defined color palette
            lwd = 2,                     # set line width
            lty = line_type[cr])                  # set line type
      
      polygon(x = c(1:time2, rev(1:time2)), 
              y = c(B_lower[a, , cr], rev(B_upper[a, , cr])), 
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
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
    dr_ytick <- seq(y1_dr, y2_dr, by_dr)              # set y axis tick marks
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
      lines(DR_medians[, cr],
            col = color[cr],                   # use pre-defined color palette
            lwd = 2,                           # set line width
            lty = (cr %% 3) + 1)}              # set line type
    
    # add a gray dotted line at y = final_target_DR
    lines(0:time2, rep(0.6, time2 + 1), col = 'gray', lty = 3)
    
    # add a legend
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(1, type = 'n', axes = F, xlab = '', ylab = '')
    legend(x = 'top', inset = 0, horiz = T,    # position
           col = color,                              # apply color palette
           lwd = 2,                                  # apply line thicknesses
           lty = line_type,                          # apply line patterns
           title = expression(bold('Control Rule')), # add legend title, labels
           c("Static \n Low M", "Static \n Correct M", "Static \n High M",
             "Transient \n Low M", "Transient \n Correct M", "Transient \n High M"),           
           seg.len = 3,                        # adjust length of lines
           cex = 1.1, 
           bty = 'n') 
    
  }
  
  # ##### Plot relative SSB + range over time after reserve implementation #####
  # 
  # # y-axis limits
  # yyy1 <- 0.5
  # yyy2 <- 6.5
  # yyy_by <- (yyy2 - yyy1)/2
  # 
  # # x-axis limits
  # x1 <- 0
  # x2 <- time2
  # x_by <- x2/2
  # 
  # for (a in 1:3) {
  # 
  #   area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
  #   title <- sprintf("Relative SSB range: %s reserve", area)
  # 
  #   # plot the relative yield
  #   plot(1, type = 'l',                          # make an empty line graph
  #        main = title,                           # title of plot
  #        ylab = 'Relative SSB',                # axis labels
  #        xlab = 'Years since marine reserve implementation',
  #        xaxt = 'n',
  #        yaxt = 'n',                             # get rid of y-axis
  #        xlim = c(x1, x2),                     # set x-axis limits
  #        ylim = c(yyy1, yyy2))
  # 
  #   # add a gray dotted line at y = 1
  #   lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3)
  # 
  #   # set specific y-axis
  #   yyytick <- seq(yyy1, yyy2, by = yyy_by)          # set yaxis tick marks
  #   axis(side = 2,                               # specify y axis
  #        at = yyytick,                            # apply tick marks
  #        labels = T,                             # apply appropriate labels
  #        las = 1)                                # set text horizontal
  # 
  #   # set specific x-axis
  #   xtick <- seq(x1, x2, by = x_by)          # set x axis tick marks
  #   axis(side = 1,                               # specify x axis
  #        at = xtick,                            # apply tick marks
  #        labels = T,                             # apply appropriate labels
  #        las = 1)                                # set text horizontal
  # 
  #   for (cr in 1:CR) {
  #       lines(SSB_medians[a, , cr],
  #             col = color[cr],                  # use pre-defined color palette
  #             lwd = 2,                     # set line width
  #             lty = line_type[cr])                  # set line type
  #       
  #       polygon(x = c(1:time2, rev(1:time2)), 
  #               y = c(SSB_lower[a, , cr], rev(SSB_upper[a, , cr])), 
  #               col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
  #   }
  #   
  #   # add a legend
  #   legend(x = c(21.25, 29), y = c(4, 6.7),   # position
  #          col = color,                           # apply viridis color palette
  #          lwd = 2,                # apply line thicknesses
  #          lty = line_type,             # apply line patterns
  #          title = expression(bold('Control Rule')),                         # add legend title and labels
  #          c("Static Low M", "Static Correct M", "Static High M", 
  #            "Transient Low M", "Transient Correct M", "Transient High M"),           seg.len = 3.5,                        # adjust length of lines
  #          cex = 1)                            # text size
  # }
  
  if (plot_individual_runs == T) {
    
    for (a in 1:2) {
      
      # y-axis limits
      y1 <- 0
      y2 <- 3
      y_by <- (y2 - y1)/4
      
      # x-axis limits
      x1 <- 0
      x2 <- time2
      x_by <- x2/4
      
      area <- ifelse(a == 1, 'far from', 'near')
      title <- sprintf("Relative yield: %s reserve", area)
      
      # plot the relative biomass
      plot(1, type = 'l',                          # make an empty line graph
           main = title,                           # title of plot
           ylab = 'Relative Yield',              # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(x1, x2),                       # set x-axis limits
           ylim = c(y1, y2), 
           cex.main = 1.75, cex.axis = 3, cex.lab = 1.5)
      
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
      
      for (k in 1:sample_size) {
        
        for (cr in c(2, 5)) {
          lines(Y_sample_runs[a, , cr, k],
                col = color[cr],                  # use pre-defined color palette
                lwd = 2,                     # set line width
                lty = line_type[cr])                  # set line type
        }
        
      }
    
    # add a legend
    legend(x = c(21.25, 28.25), y = c(1.75, 3.11),   # position
           col = color,                           # apply viridis color palette
           lwd = 2,                # apply line thicknesses
           lty = line_type,             # apply line patterns
           title = expression(bold('Control Rule')),                         # add legend title and labels
           c("Static Low M", "Static Correct M", "Static High M", 
             "Transient Low M", "Transient Correct M", "Transient High M"),           
           seg.len = 3.5,                        # adjust length of lines
           cex = 1)   
    
    }
    
  }
  
}