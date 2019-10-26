# plot_stuff <- function(yield.in, biomass.in, A, time2, CR) {
  
setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")  

library(viridis)
  
  load("../../data/sims_yield.Rda")
  load("../../data/sims_biomass.Rda")
  
  # replace NA's in sims_yield with 0s
  sims_yield[is.na(sims_yield)] <- 0
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  Y_lowerIQR <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  Y_upperIQR <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  
  B_medians <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  B_lowerIQR <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  B_upperIQR <- array(rep(NA, A*time2*CR), c(A, time2, CR))
  
  # extract data from files and plot medians + interquartile ranges
  for (i in 1:A) {
    for (j in 1:time2) {
      for (k in 1:CR) {
        
        Y_medians[i, j, k] <- median(sims_yield[i, j, k, ])
        B_medians[i, j, k] <- median(sims_biomass[i, j, k, ])
        
        Y_lowerIQR[i, j, k] <- quantile(sims_yield[i, j, k, ])[[2]]
        B_lowerIQR[i, j, k] <- quantile(sims_biomass[i, j, k, ])[[2]]
        
        Y_upperIQR[i, j, k] <- quantile(sims_yield[i, j, k, ])[[4]]
        B_upperIQR[i, j, k] <- quantile(sims_biomass[i, j, k, ])[[4]]
        
      }
    }
  } 
  
  ###### Plot relative yield + IQR over time after reserve implementation ######  
  
  ### set control rules to be compared:
  comparing <- c(2, 5)
  
  # use colorblind color palette, viridis
  color <- viridis(CR)
  
  # set plot margins to leave room for legend
  par(mar = c(5.1, 4.1, 4.1, 8.7), xpd = T)
  
  # y-axis limits
  y1 <- 0.5
  y2 <- 1.5
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
    
    for (cr in comparing) {
      lines(Y_medians[a, , cr],
            col = color[cr],                     # use pre-defined color palette
            lwd = cr%%2 + 1,                     # set line width
            lty = ceiling(cr/2)                  # set line type
      )
      
      polygon(x = c(1:time2, rev(1:time2)), 
              y = c(Y_lowerIQR[a, , cr], rev(Y_upperIQR[a, , cr])), 
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
    # add a legend
    legend(x = c(54, 62), y = c(y2 + 0.04, y2 - 0.475),   # position
           col = color,                          # apply viridis color palette
           lwd = rep(c(2, 1), 4),                # apply line thicknesses
           lty = rep(1:5, each = 2),             # apply line patterns
           title = 'CR',                         # add legend title and labels
           c("CR 1", "CR 2", "CR 3", "CR 4", "CR 5", "CR 6", "CR 7", "CR 8"),
           seg.len = 3.5,                        # adjust length of lines
           cex = 0.9)                            # text size
  }
  
  ##### Plot relative biomass + IQR over time after reserve implementation #####
  
  # y-axis limits
  yy1 <- 0.5
  yy2 <- 1.5
  yy_by <- (yy2 - yy1)/2
  
  # x-axis limits
  xx1 <- 0
  xx2 <- time2
  xx_by <- xx2/2
  
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
    
    for (cr in comparing) {                           # add one line per control rule
      lines(B_medians[a, , cr],
            col = color[cr],                     # use pre-defined color palette
            lwd = cr%%2 + 1,                     # set line width
            lty = ceiling(cr/2)                  # set line type
      )
      
      polygon(x = c(1:time2, rev(1:time2)), 
              y = c(B_lowerIQR[a, , cr], rev(B_upperIQR[a, , cr])), 
              col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
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
  
# }

