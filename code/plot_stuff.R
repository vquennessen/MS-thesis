plot_stuff <- function(filepath1, filepath2, filepath3, filepath4, 
                       A, Time1, Time2, CR, num_sims, sample_size, PD, 
                       plot_individual_runs, y_DR, species, final_DR) {
  
  # setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")  
  
  # load objects
  load(filepath1)
  load(filepath2)
  load(filepath3)
  load(filepath4)
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  # pull out sample sims
  Y_sample <- sims_yield[, , , indices]
  B_sample <- sims_biomass[, , , indices]
  SSB_sample <- sims_SSB[, , , indices]
  DR_sample <- sims_DR[, , indices]
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, A*(Time2 + 1)*CR*num_sims), c(A, Time2 + 1, CR, num_sims))
  Rel_yield <- array(rep(0, A*(Time2 + 1)*CR*num_sims), c(A, Time2 + 1, CR, num_sims))
  Rel_SSB <- array(rep(0, A*(Time2 + 1)*CR*num_sims), c(A, Time2 + 1, CR, num_sims))
  
  # total time
  TimeT <- Time1 + Time2
  
  # calculate relative arrays after reserve implementation
  for (a in 1:A) {
    for (cr in 1:CR) {
      for (sim in indices) {
        Rel_biomass[a, , cr, sim] <- sims_biomass[a, Time1:TimeT, cr, ENM, sim]/sims_biomass[a, Time1, cr, ENM, sim]
        Rel_yield[a, , cr, sim] <- sims_yield[a, Time1:TimeT, cr, ENM, sim]/sims_yield[a, Time1, cr, ENM, sim]
        Rel_SSB[a, , cr, sim] <- sims_SSB[a, Time1:TimeT, cr, ENM, sim]/sims_SSB[a, Time1, cr, ENM, sim]        
      }
    }
  }
  
  # replace NA's in sims_yield with 0s
  Y_sample_runs[is.na(Y_sample_runs)] <- 0
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- Y_lower <- Y_upper <- array(rep(NA, A*(time2 + 1)*CR), 
                                           c(A, time2 + 1, CR))
  B_medians <- B_lower <- B_upper <- array(rep(NA, A*(time2 + 1)*CR), 
                                           c(A, time2 + 1, CR))
  SSB_medians <- SSB_lower <- SSB_upper <- array(rep(NA, A*(time2 + 1)*CR), 
                                                 c(A, time2 + 1, CR))
  DR_medians <- array(rep(NA, (time2 + 1)*CR), c(time2 + 1, CR))
  
  # extract data from files and plot medians + interquartile ranges
  for (t in 1:(time2 + 1)) {
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
  
  # use red-blue color palette
  palette <- colorRampPalette(c('red', 'blue'))
  color <- palette(CR)
  
  # set line types - solid for correct M, dashed for high M, dotted for low M
  line_type <- c(2, 1, 3, 2, 1, 3)
  
  # set layout matrix for all plots
  layout_m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = T)
  w1 <- 2; w2 <- 0.35
  h1 <- 5; h2 <- 2
  
  # set main title
  main_title <- paste(species, ", Target DR = ", final_DR, sep = '')
  
  # set legend title and text and position
  legend_title <- expression(bold('Control Rule'))
  legend_text  <- c("\n Static \n Low M", "\n Static \n Correct M", 
                    "\n Static \n High M", "\n Transient \n Low M", 
                    "\n Transient \n Correct M", "\n Transient \n High M")
  position <- 'left'
  
  # plot margins
  plot1_margins <- c(0.75, 4.7, 6, 0.25)
  plot2_margins <- c(4.3, 4.7, 2, 0.25)
  plot3_margins <- c(0, 1, 0, 0)
  
  # text sizes
  mt <- 1.75    # main title for main plot
  mt2 <- 1.5    # main title for DR plot
  lab <- 1.5    # axis labels
  ax <- 1.25    # axis tick labels
  leg <- 1.25   # legend text
  
  ##### Plot relative biomass + range over time after reserve implementation #####
  
  # y-axis limits
  y1 <- 0
  y2 <- 4
  y_by <- (y2 - y1)/4
  
  # x-axis limits
  x1 <- 0
  x2 <- time2
  x_by <- x2/4
  
  # DR y-axis limits
  y1_dr <- 0
  y2_dr <- 2
  by_dr <- y2_dr/2
  
  for (a in 1:3) {
    
    # set plotting layout
    layout(mat = layout_m,
           widths = c(w1, w2),                   # Widths of the 2 columns
           heights = c(h1, h2))                  # Heights of the 2 rows
    
    area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
    sub_title <- sprintf("Relative biomass: %s reserve", area)
    
    # plot the relative biomass
    par(mar = plot1_margins)
    plot(1, type = 'l',                          # make an empty line graph
         main = paste(main_title, '\n', sub_title, sep = ''),  # title of plot
         ylab = 'Relative Biomass',              # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                       # set x-axis limits
         ylim = c(y1, y2), 
         cex.main = mt, cex.lab = lab)
    
    # add a gray dotted line at y = 1
    lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3, lwd = 2)
    
    # set specific y-axis
    ytick <- seq(y1, y2, by = y_by)              # set yaxis tick marks
    axis(side = 2,                               # specify y axis
         at = ytick,                             # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1,                                # set text horizontal
         cex.axis = ax)                          # set axis text size 
    
    for (cr in 1:CR) {
      lines(x1:x2, B_medians[a, , cr],
            col = color[cr],                     # use pre-defined color palette
            lwd = 2,                             # set line width
            lty = line_type[cr])                 # set line type
      
      polygon(x = c(x1:x2, rev(x1:x2)), 
              y = c(B_lower[a, , cr], rev(B_upper[a, , cr])), 
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
    # plot the density ratio over time
    par(mar = plot2_margins)
    plot(1, type = 'l',                          # make an empty line graph
         main = 'Density Ratios Over Time',      # title of plot
         ylab = 'Density Ratio',                 # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(0, time2),                     # set x-axis limits
         ylim = c(0, y2_dr), 
         cex.lab = lab, cex.main = mt2)
    
    # set specific y-axis
    dr_ytick <- seq(y1_dr, y2_dr, by_dr)         # set y axis tick marks
    axis(side = 2,                               # specify y axis
         at = dr_ytick,                          # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1,                                # set text horizontal
         cex.axis = ax)                          # set axis text size 
    
    # set specific x-axis
    xtick <- seq(x1, x2, by = x_by)              # set x axis tick marks
    axis(side = 1,                               # specify x axis
         at = xtick,                             # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1,                                # set text horizontal
         cex.axis = ax)                          # set axis text size 
    
    for (cr in 1:CR) {
      lines(x1:x2, DR_medians[, cr],
            col = color[cr],                     # use pre-defined color palette
            lwd = 2,                             # set line width
            lty = (cr %% 3) + 1)}                # set line type
    
    # add a gray dotted line at target_DR over time
    lines(0:time2, y_DR, col = 'gray', lty = 3, lwd = 2)
    
    # add a legend
    par(mar = plot3_margins)
    plot(1, type = 'n', axes = F, xlab = '', ylab = '')
    legend(x = position, inset = 0, horiz = F,   # position
           col = color,                          # apply color palette
           lwd = 2,                              # apply line thicknesses
           lty = line_type,                      # apply line patterns
           title = legend_title,                 # add legend title
           legend = legend_text,                 # add legend labels
           seg.len = 3,                          # adjust length of lines
           cex = leg,                            # adjust legend text size
           bty = 'n') 
    
  }

  ###### Plot relative yield + IQR over time after reserve implementation ######

  # y-axis limits
  yy1 <- 0.3
  yy2 <- 1.5
  yy_by <- (yy2 - yy1)/4

  for (a in 1:2) {

    area <- ifelse(a == 1, 'far from', 'near')
    sub_title <- sprintf("Relative yield: %s reserve", area)

    # set plotting layout
    layout(mat = layout_m,
           widths = c(w1, w2),                   # Widths of the 2 columns
           heights = c(h1, h2))                  # Heights of the 2 rows

    # plot the relative yield
    par(mar = plot1_margins)
    plot(1, type = 'l',                          # make an empty line graph
         main = paste(main_title, '\n', sub_title, sep = ''),  # title of plot
         ylab = 'Relative Yield',                # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                       # set x-axis limits
         ylim = c(yy1, yy2),
         cex.main = mt, cex.lab = lab)

    # add a gray dotted line at y = 1
    lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3, lwd = 2)

    # set specific y-axis
    yytick <- seq(yy1, yy2, by = yy_by)        # set y axis tick marks
    axis(side = 2,                             # specify y axis
         at = yytick,                          # apply tick marks
         labels = T,                           # apply appropriate labels
         las = 1,                              # set text horizontal
         cex.axis = ax)                        # set axis text size 
    
    for (cr in 1:CR) {
      lines(x1:x2, Y_medians[a, , cr],
            col = color[cr],                   # use pre-defined color palette
            lwd = 2,                           # set line width
            lty = line_type[cr])               # set line type

      polygon(x = c(x1:x2, rev(x1:x2)),
              y = c(Y_lower[a, , cr], rev(Y_upper[a, , cr])),
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }

    # plot the density ratios over time
    par(mar = plot2_margins)
    plot(1, type = 'l',                        # make an empty line graph
         main = 'Density Ratios Over Time',    # title of plot
         ylab = 'Density Ratio',               # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                           # get rid of y-axis
         xlim = c(0, time2),                   # set x-axis limits
         ylim = c(0, y2_dr),
         cex.lab = lab, cex.main = mt2)

    # set specific y-axis
    dr_ytick <- seq(y1_dr, y2_dr, by_dr)       # set y axis tick marks
    axis(side = 2,                             # specify y axis
         at = dr_ytick,                        # apply tick marks
         labels = T,                           # apply appropriate labels
         las = 1,                              # set text horizontal
         cex.axis = ax)                        # set axis text size 
    
    # set specific x-axis
    xtick <- seq(x1, x2, by = x_by)            # set x axis tick marks
    axis(side = 1,                             # specify x axis
         at = xtick,                           # apply tick marks
         labels = T,                           # apply appropriate labels
         las = 1,                              # set text horizontal
         cex.axis = ax)                        # set axis text size 
    
    for (cr in 1:CR) {
      lines(x1:x2, DR_medians[, cr],
            col = color[cr],                   # use pre-defined color palette
            lwd = 2,                           # set line width
            lty = (cr %% 3) + 1)}              # set line type

    # add a gray dotted line at target_DR over time
    lines(0:time2, y_DR, col = 'gray', lty = 3, lwd = 2)
    
    # add a legend
    par(mar = plot3_margins)
    plot(1, type = 'n', axes = F, xlab = '', ylab = '')
    legend(x = position, inset = 0, horiz = F,   # position
           col = color,                          # apply color palette
           lwd = 2,                              # apply line thicknesses
           lty = line_type,                      # apply line patterns
           title = legend_title,                 # add legend title
           legend = legend_text,                 # add legend labels
           seg.len = 3,                          # adjust length of lines
           cex = leg,                            # adjust legend text size
           bty = 'n') 
  }

  ##### Plot relative SSB + range over time after reserve implementation #####
  
  # y-axis limits
  yyy1 <- 0.5
  yyy2 <- 6.5
  yyy_by <- (yyy2 - yyy1)/4
  
  for (a in 1:3) {
  
    area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
    sub_title <- sprintf("Relative SSB: %s reserve", area)
    
    # set plotting layout
    layout(mat = layout_m,
           widths = c(w1, w2),                   # Widths of the 2 columns
           heights = c(h1, h2))                  # Heights of the 2 rows
  
    # plot the relative SSB
    par(mar = plot1_margins)
    plot(1, type = 'l',                          # make an empty line graph
         main = paste(main_title, '\n', sub_title, sep = ''),  # title of plot
         ylab = 'Relative SSB',                  # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                       # set x-axis limits
         ylim = c(yyy1, yyy2), 
         cex.main = mt, cex.lab = lab)
    
    # add a gray dotted line at y = 1
    lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3, lwd = 2)
  
    # set specific y-axis
    yyytick <- seq(yyy1, yyy2, by = yyy_by)      # set yaxis tick marks
    axis(side = 2,                               # specify y axis
         at = yyytick,                           # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1,                                # set text horizontal
         cex.axis = ax)                          # set axis text size             
  
    for (cr in 1:CR) {
        lines(x1:x2, SSB_medians[a, , cr],
              col = color[cr],                   # use pre-defined color palette
              lwd = 2,                           # set line width
              lty = line_type[cr])               # set line type
  
        polygon(x = c(x1:x2, rev(x1:x2)),
                y = c(SSB_lower[a, , cr], rev(SSB_upper[a, , cr])),
                col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
   }
  
    # plot the density ratios over time
    par(mar = plot2_margins)
    plot(1, type = 'l',                        # make an empty line graph
         main = 'Density Ratios Over Time',    # title of plot
         ylab = 'Density Ratio',               # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                           # get rid of y-axis
         xlim = c(0, time2),                   # set x-axis limits
         ylim = c(0, y2_dr),
         cex.lab = lab, cex.main = mt2)
    
    # set specific y-axis
    dr_ytick <- seq(y1_dr, y2_dr, by_dr)       # set y axis tick marks
    axis(side = 2,                             # specify y axis
         at = dr_ytick,                        # apply tick marks
         labels = T,                           # apply appropriate labels
         las = 1,                              # set text horizontal
         cex.axis = ax)                        # set axis text size
    
    # set specific x-axis
    xtick <- seq(x1, x2, by = x_by)            # set x axis tick marks
    axis(side = 1,                             # specify x axis
         at = xtick,                           # apply tick marks
         labels = T,                           # apply appropriate labels
         las = 1,                              # set text horizontal
         cex.axis = ax)                        # set axis text size             

    
    for (cr in 1:CR) {
      lines(x1:x2, DR_medians[, cr],
            col = color[cr],                   # use pre-defined color palette
            lwd = 2,                           # set line width
            lty = (cr %% 3) + 1)}              # set line type
    
    # add a gray dotted line at target_DR over time
    lines(0:time2, y_DR, col = 'gray', lty = 3, lwd = 2)
    
    # add a legend
    par(mar = plot3_margins)
    plot(1, type = 'n', axes = F, xlab = '', ylab = '')
    legend(x = position, inset = 0, horiz = F,   # position
           col = color,                          # apply color palette
           lwd = 2,                              # apply line thicknesses
           lty = line_type,                      # apply line patterns
           title = legend_title,                 # add legend title
           legend = legend_text,                 # add legend labels
           seg.len = 3,                          # adjust length of lines
           cex = leg,                            # adjust legend text size
           bty = 'n') 
  }
  
  ##### Plot individual runs ###################################################

  if (plot_individual_runs == T) {

    ##### Yield #####

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
           ylab = 'Relative Yield',                # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(x1, x2),                       # set x-axis limits
           ylim = c(y1, y2),
           cex.main = 1.75, cex.axis = 3, cex.lab = 1.5)

      # add a gray dotted line at y = 1
      lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3, lwd = 2)

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
          lines(x1:x2, Y_sample_runs[a, , cr, k],
                col = color[cr],                  # use pre-defined color palette
                lwd = 2,                     # set line width
                lty = line_type[cr])                  # set line type
        }

      }

      # add a legend
      legend(x = 'right', inset = 0, horiz = F,   # position
             col = color,                          # apply color palette
             lwd = 2,                              # apply line thicknesses
             lty = line_type,                      # apply line patterns
             title = legend_title,                 # add legend title
             legend = legend_text,                 # add legend labels
             seg.len = 3,                          # adjust length of lines
             cex = 1.1,                            # adjust legend text size
             bty = 'n') 
    }

    
    ##### Biomass #####
    
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
      title <- sprintf("Relative biomass: %s reserve", area)
      
      # plot the relative biomass
      plot(1, type = 'l',                          # make an empty line graph
           main = title,                           # title of plot
           ylab = 'Relative Biomass',              # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(x1, x2),                       # set x-axis limits
           ylim = c(y1, y2),
           cex.main = 1.75, cex.axis = 3, cex.lab = 1.5)
      
      # add a gray dotted line at y = 1
      lines(0:time2, rep(1, time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
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
          lines(x1:x2, B_sample_runs[a, , cr, k],
                col = color[cr],                  # use pre-defined color palette
                lwd = 2,                     # set line width
                lty = line_type[cr])                  # set line type
        }
        
      }
      
      # add a legend
      legend(x = 'right', inset = 0, horiz = F,   # position
             col = color,                          # apply color palette
             lwd = 2,                              # apply line thicknesses
             lty = line_type,                      # apply line patterns
             title = legend_title,                 # add legend title
             legend = legend_text,                 # add legend labels
             seg.len = 3,                          # adjust length of lines
             cex = 1.1,                            # adjust legend text size
             bty = 'n') 
    }
    
  }

}