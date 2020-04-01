  ##### set up #####
  # setwd("C:/Users/Vic/Documents/Projects/DensityRatio/code/R")  
  library(densityratio)
  
  # set variables
  Species = 'CAB_OR_2019'
  num_sims = 2

  # # load objects
  # load(paste('~/Projects/MS-thesis/data/', Species, '/', num_sims, 
  #            '_N.Rda', sep = ''))
  # load(paste('~/Projects/MS-thesis/data/', Species, '/', num_sims, 
  #            '_biomass.Rda', sep = ''))
  # load(paste('~/Projects/MS-thesis/data/', Species, '/', num_sims, 
  #            '_SSB.Rda', sep = ''))
  # load(paste('~/Projects/MS-thesis/data/', Species, '/', num_sims, 
  #            '_yield.Rda', sep = ''))
  # load(paste('~/Projects/MS-thesis/data/', Species, '/', num_sims, 
  #            '_effort.Rda', sep = ''))
  # load(paste('~/Projects/MS-thesis/data/', Species, '/', num_sims, 
  #            '_DR.Rda', sep = ''))
  # 
  # set other variables (should never really change for now)
  A = 5
  MPA = 3
  Time1 = 50
  Time2 = 20
  Final_DRs = c(0.2, 0.4, 0.6, 0.8, 1)
  Control_rules = c(1:6)
  
  # dimensions
  TimeT <- Time1 + Time2
  CR <- length(Control_rules)
  FDR <- length(Final_DRs)
  Rec_age <- parameters(Species)[[3]]
  Max_age <- parameters(Species)[[1]]
  ages <- Rec_age:Max_age
  n <- length(ages)
  
  # plotting variables
  sample_size = num_sims
  PD = 0.25
  plot_individual_runs = FALSE
  Nat_mortality <- parameters(Species)[[2]]
  Error = 0.05
  ENM = 2
  
  # sample from simulations
  indices <- sample(1:num_sims, sample_size, replace = FALSE)
  
  ##### process yield, biomass, SSB, effort, and DR   #####
  
  # pull out sample sims
  Y_sample   <- sims_yield[, , , , indices]
  B_sample   <- sims_biomass[, , , , indices]
  SSB_sample <- sims_SSB[, , , , indices]
  E_sample   <- sims_effort[, , , indices]
  DR_sample  <- sims_DR[, , , indices]
  A_sample   <- sims_abundance[, , , , indices]
  
  # initialize relative arrays
  Rel_biomass <- array(rep(0, MPA*(Time2 + 1)*CR*FDR*num_sims), 
                       c(MPA, Time2 + 1, CR, FDR, num_sims))
  Rel_yield <- array(rep(0, MPA*(Time2 + 1)*CR*FDR*num_sims), 
                     c(MPA, Time2 + 1, CR, FDR, num_sims))
  Rel_SSB <- array(rep(0, MPA*(Time2 + 1)*CR*FDR*num_sims), 
                   c(MPA, Time2 + 1, CR, FDR, num_sims))
  Rel_effort <- array(rep(0, (Time2 + 1)*CR*FDR*num_sims), 
                      c(Time2 + 1, CR, FDR, num_sims))
  Rel_abundance <- array(rep(0, MPA*(Time2 + 1)*CR*FDR*num_sims), 
                   c(MPA, Time2 + 1, CR, FDR, num_sims))
  
  # calculate relative arrays after reserve implementation
  for (cr in 1:CR) {
    for (fdr in 1:FDR) {
      for (sim in indices) {
        Rel_effort[ , cr, fdr, sim] <- E_sample[Time1:TimeT, cr, fdr, sim] /
          E_sample[Time1, cr, fdr, sim]        
        for (a in 1:(MPA - 1)) {
          Rel_biomass[a, , cr, fdr, sim] <- B_sample[a, , cr, fdr, sim] / 
            B_sample[a, 1, cr, fdr, sim]
          Rel_SSB[a, , cr, fdr, sim] <- SSB_sample[a, , cr, fdr, sim] / 
            SSB_sample[a, 1, cr, fdr, sim]  
          Rel_abundance[a, , cr, fdr, sim] <- SSB_sample[a, , cr, fdr, sim] / 
            SSB_sample[a, 1, cr, fdr, sim] 
          Rel_yield[a, , cr, fdr, sim] <- Y_sample[a, , cr, fdr, sim] / 
            Y_sample[a, 1, cr, fdr, sim] 
          for (a in MPA) {
            Rel_biomass[a, , cr, fdr, sim] <- B_sample[a, , cr, fdr, sim] / 
              B_sample[a, 1, cr, fdr, sim]
            Rel_SSB[a, , cr, fdr, sim] <- SSB_sample[a, , cr, fdr, sim] / 
              SSB_sample[a, 1, cr, fdr, sim]  
            Rel_abundance[a, , cr, fdr, sim] <- SSB_sample[a, , cr, fdr, sim] / 
              SSB_sample[a, 1, cr, fdr, sim] 
          }
        }
      }
    }
  }
  
  # initialize median, lowerIQR, and upperIQR arrays
  Y_medians <- Y_lower <- Y_upper <- array(rep(NA, (MPA - 1)*(Time2 + 1)*CR*FDR), 
                                           c(MPA - 1, Time2 + 1, CR, FDR))
  B_medians <- B_lower <- B_upper <- array(rep(NA, MPA*(Time2 + 1)*CR*FDR), 
                                           c(MPA, Time2 + 1, CR, FDR))
  SSB_medians <- SSB_lower <- SSB_upper <- array(rep(NA, MPA*(Time2 + 1)*CR*FDR), 
                                                 c(MPA, Time2 + 1, CR, FDR))
  A_medians <- A_lower <- A_upper <- array(rep(NA, MPA*(Time2 + 1)*CR*FDR), 
                                                 c(MPA, Time2 + 1, CR, FDR))
  DR_medians <- array(rep(NA, (Time2 + 1)*CR*FDR), c(Time2 + 1, CR, FDR))
  E_medians <- E_lower <- E_upper <- array(rep(NA, (Time2 + 1)*CR*FDR), 
                                           c(Time2 + 1, CR, FDR))
  # calculate medians, upper limits, and lower limits  
  for (t in 1:(Time2 + 1)) {
    for (cr in 1:CR) {
      for (fdr in 1:FDR) {        
        
        DR_medians[t, cr, fdr] <- median(DR_sample[t + Time1 - 1, cr, fdr, ])
        
        E_medians[t, cr, fdr] <- median(Rel_effort[t, cr, fdr, ])
        E_lower[t, cr, fdr] <- quantile(Rel_effort[t, cr, fdr, ], 0.5 - PD)
        E_upper[t, cr, fdr] <- quantile(Rel_effort[t, cr, fdr, ], 0.5 + PD)        
        
        for (a in 1:(MPA - 1)) {
          B_medians[a, t, cr, fdr] <- median(Rel_biomass[a, t, cr, fdr, ])
          SSB_medians[a, t, cr, fdr] <- median(Rel_SSB[a, t, cr, fdr, ])
          A_medians[a, t, cr, fdr] <- median(Rel_abundance[a, t, cr, fdr, ])
          Y_medians[a, t, cr, fdr] <- median(Rel_yield[a, t, cr, fdr, ])
          
          B_lower[a, t, cr, fdr] <- quantile(Rel_biomass[a, t, cr, fdr, ], 0.5 - PD)
          SSB_lower[a, t, cr, fdr] <- quantile(Rel_SSB[a, t, cr, fdr, ], 0.5 - PD)
          A_lower[a, t, cr, fdr] <- quantile(Rel_abundance[a, t, cr, fdr, ], 0.5 - PD)
          Y_lower[a, t, cr, fdr] <- quantile(Rel_yield[a, t, cr, fdr, ], 0.5 - PD)
          
          B_upper[a, t, cr, fdr] <- quantile(Rel_biomass[a, t, cr, fdr, ], 0.5 + PD)
          SSB_upper[a, t, cr, fdr] <- quantile(Rel_SSB[a, t, cr, fdr, ], 0.5 + PD)
          A_upper[a, t, cr, fdr] <- quantile(Rel_abundance[a, t, cr, fdr, ], 0.5 + PD)
          Y_upper[a, t, cr, fdr] <- quantile(Rel_yield[a, t, cr, fdr, ], 0.5 + PD)
          
        }
        
        for (a in MPA) {
          
          B_medians[a, t, cr, fdr] <- median(Rel_biomass[a, t, cr, fdr, ])
          SSB_medians[a, t, cr, fdr] <- median(Rel_SSB[a, t, cr, fdr, ])
          A_medians[a, t, cr, fdr] <- median(Rel_abundance[a, t, cr, fdr, ])
          
          B_lower[a, t, cr, fdr] <- quantile(Rel_biomass[a, t, cr, fdr, ], 0.5 - PD)
          SSB_lower[a, t, cr, fdr] <- quantile(Rel_SSB[a, t, cr, fdr, ], 0.5 - PD)
          A_lower[a, t, cr, fdr] <- quantile(Rel_abundance[a, t, cr, fdr, ], 0.5 - PD)
          
          B_upper[a, t, cr, fdr] <- quantile(Rel_biomass[a, t, cr, fdr, ], 0.5 + PD)
          SSB_upper[a, t, cr, fdr] <- quantile(Rel_SSB[a, t, cr, fdr, ], 0.5 + PD)
          A_upper[a, t, cr, fdr] <- quantile(Rel_abundance[a, t, cr, fdr, ], 0.5 + PD)
          
        }
      }
    }
  } 
  
  #####
  ##### plotting parameters (biomass) #####
  
  # use red-blue color palette
  palette <- colorRampPalette(c('red', 'blue'))
  color <- palette(FDR)
  
  # set line types - solid for correct M, dashed for high M, dotted for low M
  line_type <- c(1:FDR)
  
  # set layout matrix for all plots
  layout_m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = T)
  w1 <- 2; w2 <- 0.5
  h1 <- 3; h2 <- 2
  
  # set legend title and text and position
  legend_title <- expression(bold('Final DR'))
  legend_text  <- c(Final_DRs)
  position <- 'left'
  
  # plot margins
  plot1_margins <- c(0.75, 5, 6, 0.25)
  plot2_margins <- c(4.3, 5, 2, 0.25)
  plot3_margins <- c(0, 1, 0, 0)
  
  # text sizes
  mt <- 1.5    # main title for main plot
  mt2 <- 1.5    # main title for DR plot
  lab <- 1.5    # axis labels
  ax <- 1.25    # axis tick labels
  leg <- 1.25   # legend text
  
  # y-axis limits
  y1 <- 0
  y2 <- 2
  y_by <- (y2 - y1)/4
  
  # x-axis limits
  x1 <- 0
  x2 <- Time2
  x_by <- x2/4
  
  # DR y-axis limits
  y1_dr <- 0.5
  y2_dr <- 1.5
  by_dr <- (y2_dr - y1_dr)/2
  
  ##### Plot relative biomass and DR across FDRs #####
  
  # set plotting layout
  layout(mat = layout_m,
         widths = c(w1, w2),                   # Widths of the 2 columns
         heights = c(h1, h2))                  # Heights of the 2 rows
  
  for (cr in c(2, 5)){
    
    # set title
    cr_label <- ifelse(cr == 2, 'Static', 'Transient')
    main_title <- paste(Species, ' - ', cr_label, ' CR', sep = '')
    
    for (a in 1:MPA) {
      
      area <- ifelse(a == 1, 'far from', ifelse(a == 2, 'near', 'in'))
      sub_title <- sprintf("Relative Biomass: %s reserve", area)
      
      # plot the relative effort
      par(mar = plot1_margins)
      plot(1, type = 'l',                          # make an empty line graph
           main = paste(main_title, '\n', sub_title, sep = ''),  # title of plot
           ylab = 'Relative Biomass',              # axis labels
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(x1, x2),                       # set x-axis limits
           ylim = c(y1, y2), 
           cex.main = mt, cex.lab = lab)
      
      # add a gray dotted line at y = 1
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
      # set specific y-axis
      ytick <- seq(y1, y2, by = y_by)              # set yaxis tick marks
      axis(side = 2,                               # specify y axis
           at = ytick,                             # apply tick marks
           labels = T,                             # apply appropriate labels
           las = 1,                                # set text horizontal
           cex.axis = ax)                          # set axis text size 
      
      for (fdr in 1:FDR) {
        lines(x1:x2, B_medians[a, , cr, fdr],
              col = color[fdr],                     # use pre-defined color palette
              lwd = 2,                             # set line width
              lty = line_type[fdr])                 # set line type
        
        # polygon(x = c(x1:x2, rev(x1:x2)),
        #         y = c(B_lower[a, , cr], rev(B_upper[a, , cr])),
        #         col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
      }
      
      # plot the density ratio over time
      par(mar = plot2_margins)
      plot(1, type = 'l',                          # make an empty line graph
           main = 'Density Ratios Over Time',      # title of plot
           ylab = 'Density Ratio',                 # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(0, Time2),                     # set x-axis limits
           ylim = c(y1_dr, y2_dr), 
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
      
      for (fdr in 1:FDR) {
        lines(x1:x2, DR_medians[, cr, fdr],
              col = color[fdr],                     # use pre-defined color palette
              lwd = 2,                             # set line width
              lty = line_type[fdr])
      }                # set line type
      
      # add a gray dotted line at y = 0.6, 0.8, 1
      lines(0:Time2, rep(0.6, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      lines(0:Time2, rep(0.8, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
      # add a legend
      par(mar = plot3_margins)
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F,   # position
             col = color,                          # apply color palette
             lwd = 2,                              # apply line thicknesses
             lty = line_type,                      # apply line patterns
             title = legend_title,                 # add legend title
             legend = legend_text,                 # add legend labels
             seg.len = 2.5,                          # adjust length of lines
             cex = leg,                            # adjust legend text size
             bty = 'n')
      
    }
  }
  
  
  
  
  #
  #####
  ##### plotting parameters (yield) #####
  
  # use red-blue color palette
  palette <- colorRampPalette(c('red', 'blue'))
  color <- palette(FDR)
  
  # set line types - solid for correct M, dashed for high M, dotted for low M
  line_type <- c(1:FDR)
  
  # set layout matrix for all plots
  layout_m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = T)
  w1 <- 2; w2 <- 0.5
  h1 <- 3; h2 <- 2
  
  # set legend title and text and position
  legend_title <- expression(bold('Final DR'))
  legend_text  <- c(Final_DRs)
  position <- 'left'
  
  # plot margins
  plot1_margins <- c(0.75, 5, 6, 0.25)
  plot2_margins <- c(4.3, 5, 2, 0.25)
  plot3_margins <- c(0, 1, 0, 0)
  
  # text sizes
  mt <- 1.5    # main title for main plot
  mt2 <- 1.5    # main title for DR plot
  lab <- 1.5    # axis labels
  ax <- 1.25    # axis tick labels
  leg <- 1.25   # legend text
  
  # y-axis limits
  y1 <- 0
  y2 <- 2
  y_by <- (y2 - y1)/4
  
  # x-axis limits
  x1 <- 0
  x2 <- Time2
  x_by <- x2/4
  
  # DR y-axis limits
  y1_dr <- 0.5
  y2_dr <- 1.5
  by_dr <- (y2_dr - y1_dr)/2
  
  ##### Plot relative yield and DR across FDRs #####
  
  # set plotting layout
  layout(mat = layout_m,
         widths = c(w1, w2),                   # Widths of the 2 columns
         heights = c(h1, h2))                  # Heights of the 2 rows
  
  for (cr in c(2, 5)){
    # set title
    cr_label <- ifelse(cr == 2, 'Static', 'Transient')
    main_title <- paste(Species, ' - ', cr_label, ' CR', sep = '')
    
    for (a in 1:(MPA - 1)) {
      
      area <- ifelse(a == 1, 'far from', 'near')
      sub_title <- sprintf("Relative Yield: %s reserve", area)
      
      # plot the relative effort
      par(mar = plot1_margins)
      plot(1, type = 'l',                          # make an empty line graph
           main = paste(main_title, '\n', sub_title, sep = ''),  # title of plot
           ylab = 'Relative Yield',              # axis labels
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(x1, x2),                       # set x-axis limits
           ylim = c(y1, y2), 
           cex.main = mt, cex.lab = lab)
      
      # add a gray dotted line at y = 1
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
      # set specific y-axis
      ytick <- seq(y1, y2, by = y_by)              # set yaxis tick marks
      axis(side = 2,                               # specify y axis
           at = ytick,                             # apply tick marks
           labels = T,                             # apply appropriate labels
           las = 1,                                # set text horizontal
           cex.axis = ax)                          # set axis text size 
      
      for (fdr in 1:FDR) {
        lines(x1:x2, Y_medians[a, , cr, fdr],
              col = color[fdr],                     # use pre-defined color palette
              lwd = 2,                             # set line width
              lty = line_type[fdr])                 # set line type
        
        # polygon(x = c(x1:x2, rev(x1:x2)),
        #         y = c(B_lower[a, , cr], rev(B_upper[a, , cr])),
        #         col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
      }
      
      # plot the density ratio over time
      par(mar = plot2_margins)
      plot(1, type = 'l',                          # make an empty line graph
           main = 'Density Ratios Over Time',      # title of plot
           ylab = 'Density Ratio',                 # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(0, Time2),                     # set x-axis limits
           ylim = c(y1_dr, y2_dr), 
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
      
      for (fdr in 1:FDR) {
        lines(x1:x2, DR_medians[, cr, fdr],
              col = color[fdr],                     # use pre-defined color palette
              lwd = 2,                             # set line width
              lty = line_type[fdr])
      }                # set line type
      
      # add a gray dotted line at y = 0.6, 0.8, 1
      lines(0:Time2, rep(0.6, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      lines(0:Time2, rep(0.8, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
      # add a legend
      par(mar = plot3_margins)
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F,   # position
             col = color,                          # apply color palette
             lwd = 2,                              # apply line thicknesses
             lty = line_type,                      # apply line patterns
             title = legend_title,                 # add legend title
             legend = legend_text,                 # add legend labels
             seg.len = 2.5,                          # adjust length of lines
             cex = leg,                            # adjust legend text size
             bty = 'n')
      
    }
  }
  
  
  #####
  ##### plotting parameters (SSB) #####
  
  # use red-blue color palette
  palette <- colorRampPalette(c('red', 'blue'))
  color <- palette(FDR)
  
  # set line types - solid for correct M, dashed for high M, dotted for low M
  line_type <- c(1:FDR)
  
  # set layout matrix for all plots
  layout_m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = T)
  w1 <- 2; w2 <- 0.5
  h1 <- 3; h2 <- 2
  
  # set legend title and text and position
  legend_title <- expression(bold('Final DR'))
  legend_text  <- c(Final_DRs)
  position <- 'left'
  
  # plot margins
  plot1_margins <- c(0.75, 5, 6, 0.25)
  plot2_margins <- c(4.3, 5, 2, 0.25)
  plot3_margins <- c(0, 1, 0, 0)
  
  # text sizes
  mt <- 1.5    # main title for main plot
  mt2 <- 1.5    # main title for DR plot
  lab <- 1.5    # axis labels
  ax <- 1.25    # axis tick labels
  leg <- 1.25   # legend text
  
  # y-axis limits
  y1 <- 0
  y2 <- 2
  y_by <- (y2 - y1)/4
  
  # x-axis limits
  x1 <- 0
  x2 <- Time2
  x_by <- x2/4
  
  ##### Plot relative SSB and DR across FDRs #####
  
  # set plotting layout
  layout(mat = layout_m,
         widths = c(w1, w2),                   # Widths of the 2 columns
         heights = c(h1, h2))                  # Heights of the 2 rows
  
  for (cr in c(2, 5)){
    # set title
    cr_label <- ifelse(cr == 2, 'Static', 'Transient')
    main_title <- paste(Species, ' - ', cr_label, ' CR', sep = '')
    
    for (a in 1:MPA) {
      
      area <- ifelse(a == 1, 'far from', ifelse(a == 2, 'near', 'in'))
      sub_title <- sprintf("Relative SSB: %s reserve", area)
      
      # plot the relative effort
      par(mar = plot1_margins)
      plot(1, type = 'l',                          # make an empty line graph
           main = paste(main_title, '\n', sub_title, sep = ''),  # title of plot
           ylab = 'Relative SSB',              # axis labels
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(x1, x2),                       # set x-axis limits
           ylim = c(y1, y2), 
           cex.main = mt, cex.lab = lab)
      
      # add a gray dotted line at y = 1
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
      # set specific y-axis
      ytick <- seq(y1, y2, by = y_by)              # set yaxis tick marks
      axis(side = 2,                               # specify y axis
           at = ytick,                             # apply tick marks
           labels = T,                             # apply appropriate labels
           las = 1,                                # set text horizontal
           cex.axis = ax)                          # set axis text size 
      
      for (fdr in 1:FDR) {
        lines(x1:x2, SSB_medians[a, , cr, fdr],
              col = color[fdr],                     # use pre-defined color palette
              lwd = 2,                             # set line width
              lty = line_type[fdr])                 # set line type
        
        # polygon(x = c(x1:x2, rev(x1:x2)),
        #         y = c(B_lower[a, , cr], rev(B_upper[a, , cr])),
        #         col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
      }
      
      # plot the density ratio over time
      par(mar = plot2_margins)
      plot(1, type = 'l',                          # make an empty line graph
           main = 'Density Ratios Over Time',      # title of plot
           ylab = 'Density Ratio',                 # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(0, Time2),                     # set x-axis limits
           ylim = c(y1_dr, y2_dr), 
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
      
      for (fdr in 1:FDR) {
        lines(x1:x2, DR_medians[, cr, fdr],
              col = color[fdr],                     # use pre-defined color palette
              lwd = 2,                             # set line width
              lty = line_type[fdr])
      }                # set line type
      
      # add a gray dotted line at y = 0.6, 0.8, 1
      lines(0:Time2, rep(0.6, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      lines(0:Time2, rep(0.8, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
      # add a legend
      par(mar = plot3_margins)
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F,   # position
             col = color,                          # apply color palette
             lwd = 2,                              # apply line thicknesses
             lty = line_type,                      # apply line patterns
             title = legend_title,                 # add legend title
             legend = legend_text,                 # add legend labels
             seg.len = 2.5,                          # adjust length of lines
             cex = leg,                            # adjust legend text size
             bty = 'n')
      
    }
  }
  
  
  #
  #####
  ##### plotting parameters (abundance) #####
  
  # use red-blue color palette
  palette <- colorRampPalette(c('red', 'blue'))
  color <- palette(FDR)
  
  # set line types - solid for correct M, dashed for high M, dotted for low M
  line_type <- c(1:FDR)
  
  # set layout matrix for all plots
  layout_m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = T)
  w1 <- 2; w2 <- 0.5
  h1 <- 3; h2 <- 2
  
  # set legend title and text and position
  legend_title <- expression(bold('Final DR'))
  legend_text  <- c(Final_DRs)
  position <- 'left'
  
  # plot margins
  plot1_margins <- c(0.75, 5, 6, 0.25)
  plot2_margins <- c(4.3, 5, 2, 0.25)
  plot3_margins <- c(0, 1, 0, 0)
  
  # text sizes
  mt <- 1.5    # main title for main plot
  mt2 <- 1.5    # main title for DR plot
  lab <- 1.5    # axis labels
  ax <- 1.25    # axis tick labels
  leg <- 1.25   # legend text
  
  # y-axis limits
  y1 <- 0
  y2 <- 2
  y_by <- (y2 - y1)/4
  
  # x-axis limits
  x1 <- 0
  x2 <- Time2
  x_by <- x2/4
  
  ##### Plot relative abundance and DR across FDRs #####
  
  # set plotting layout
  layout(mat = layout_m,
         widths = c(w1, w2),                   # Widths of the 2 columns
         heights = c(h1, h2))                  # Heights of the 2 rows
  
  for (cr in c(2, 5)){
    # set title
    cr_label <- ifelse(cr == 2, 'Static', 'Transient')
    main_title <- paste(Species, ' - ', cr_label, ' CR', sep = '')
    
    for (a in 1:MPA) {
      
      area <- ifelse(a == 1, 'far from', ifelse(a == 2, 'near', 'in'))
      sub_title <- sprintf("Relative Abundance: %s reserve", area)
      
      # plot the relative effort
      par(mar = plot1_margins)
      plot(1, type = 'l',                          # make an empty line graph
           main = paste(main_title, '\n', sub_title, sep = ''),  # title of plot
           ylab = 'Relative Abundance',              # axis labels
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(x1, x2),                       # set x-axis limits
           ylim = c(y1, y2), 
           cex.main = mt, cex.lab = lab)
      
      # add a gray dotted line at y = 1
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
      # set specific y-axis
      ytick <- seq(y1, y2, by = y_by)              # set yaxis tick marks
      axis(side = 2,                               # specify y axis
           at = ytick,                             # apply tick marks
           labels = T,                             # apply appropriate labels
           las = 1,                                # set text horizontal
           cex.axis = ax)                          # set axis text size 
      
      for (fdr in 1:FDR) {
        lines(x1:x2, A_medians[a, , cr, fdr],
              col = color[fdr],                     # use pre-defined color palette
              lwd = 2,                             # set line width
              lty = line_type[fdr])                 # set line type
        
        # polygon(x = c(x1:x2, rev(x1:x2)),
        #         y = c(B_lower[a, , cr], rev(B_upper[a, , cr])),
        #         col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
      }
      
      # plot the density ratio over time
      par(mar = plot2_margins)
      plot(1, type = 'l',                          # make an empty line graph
           main = 'Density Ratios Over Time',      # title of plot
           ylab = 'Density Ratio',                 # axis labels
           xlab = 'Years since marine reserve implementation',
           xaxt = 'n',
           yaxt = 'n',                             # get rid of y-axis
           xlim = c(0, Time2),                     # set x-axis limits
           ylim = c(y1_dr, y2_dr), 
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
      
      for (fdr in 1:FDR) {
        lines(x1:x2, DR_medians[, cr, fdr],
              col = color[fdr],                     # use pre-defined color palette
              lwd = 2,                             # set line width
              lty = line_type[fdr])
      }                # set line type
      
      # add a gray dotted line at y = 0.6, 0.8, 1
      lines(0:Time2, rep(0.6, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      lines(0:Time2, rep(0.8, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
      # add a legend
      par(mar = plot3_margins)
      plot(1, type = 'n', axes = F, xlab = '', ylab = '')
      legend(x = position, inset = 0, horiz = F,   # position
             col = color,                          # apply color palette
             lwd = 2,                              # apply line thicknesses
             lty = line_type,                      # apply line patterns
             title = legend_title,                 # add legend title
             legend = legend_text,                 # add legend labels
             seg.len = 2.5,                          # adjust length of lines
             cex = leg,                            # adjust legend text size
             bty = 'n')
      
    }
  }
  
  
  #
  #####
  ##### plotting parameters (effort) #####
  
  # use red-blue color palette
  palette <- colorRampPalette(c('red', 'blue'))
  color <- palette(FDR)
  
  # set line types - solid for correct M, dashed for high M, dotted for low M
  line_type <- c(1:FDR)
  
  # set layout matrix for all plots
  layout_m <- matrix(c(1, 3, 2, 3), nrow = 2, ncol = 2, byrow = T)
  w1 <- 2; w2 <- 0.5
  h1 <- 3; h2 <- 2
  
  # set legend title and text and position
  legend_title <- expression(bold('Final DR'))
  legend_text  <- c(Final_DRs)
  position <- 'left'
  
  # plot margins
  plot1_margins <- c(0.75, 5, 6, 0.25)
  plot2_margins <- c(4.3, 5, 2, 0.25)
  plot3_margins <- c(0, 1, 0, 0)
  
  # text sizes
  mt <- 1.5    # main title for main plot
  mt2 <- 1.5    # main title for DR plot
  lab <- 1.5    # axis labels
  ax <- 1.25    # axis tick labels
  leg <- 1.25   # legend text
  
  # y-axis limits
  y1 <- 0
  y2 <- 6
  y_by <- (y2 - y1)/4
  
  # x-axis limits
  x1 <- 0
  x2 <- Time2
  x_by <- x2/4
  
  ##### Plot relative effort and DR across FDRs #####
  
  # set plotting layout
  layout(mat = layout_m,
         widths = c(w1, w2),                   # Widths of the 2 columns
         heights = c(h1, h2))                  # Heights of the 2 rows
  
  for (cr in c(2, 5)){
    
    # set title
    cr_label <- ifelse(cr == 2, 'Static', 'Transient')
    main_title <- paste(Species, ' - ', cr_label, ' CR', sep = '')
    sub_title <- sprintf('Relative Effort')
    
    # plot the relative effort
    par(mar = plot1_margins)
    plot(1, type = 'l',                          # make an empty line graph
         main = paste(main_title, '\n', sub_title, sep = ''),  # title of plot
         ylab = 'Relative Effort',              # axis labels
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(x1, x2),                       # set x-axis limits
         ylim = c(y1, y2), 
         cex.main = mt, cex.lab = lab)
    
    # add a gray dotted line at y = 1
    lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
    
    # set specific y-axis
    ytick <- seq(y1, y2, by = y_by)              # set yaxis tick marks
    axis(side = 2,                               # specify y axis
         at = ytick,                             # apply tick marks
         labels = T,                             # apply appropriate labels
         las = 1,                                # set text horizontal
         cex.axis = ax)                          # set axis text size 
    
    for (fdr in 1:FDR) {
      lines(x1:x2, E_medians[, cr, fdr],
            col = color[fdr],                     # use pre-defined color palette
            lwd = 2,                             # set line width
            lty = line_type[fdr])                 # set line type
      
      # polygon(x = c(x1:x2, rev(x1:x2)),
      #         y = c(B_lower[a, , cr], rev(B_upper[a, , cr])),
      #         col = adjustcolor(color[cr], alpha.f = 0.10), border = NA)
    }
    
    # plot the density ratio over time
    par(mar = plot2_margins)
    plot(1, type = 'l',                          # make an empty line graph
         main = 'Density Ratios Over Time',      # title of plot
         ylab = 'Density Ratio',                 # axis labels
         xlab = 'Years since marine reserve implementation',
         xaxt = 'n',
         yaxt = 'n',                             # get rid of y-axis
         xlim = c(0, Time2),                     # set x-axis limits
         ylim = c(y1_dr, y2_dr), 
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
    
    for (fdr in 1:FDR) {
      lines(x1:x2, DR_medians[, cr, fdr],
            col = color[fdr],                     # use pre-defined color palette
            lwd = 2,                             # set line width
            lty = line_type[fdr])
    }                # set line type
    
    # add a gray dotted line at y = 0.6, 0.8, 1
    lines(0:Time2, rep(0.6, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
    lines(0:Time2, rep(0.8, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
    lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
    
    # add a legend
    par(mar = plot3_margins)
    plot(1, type = 'n', axes = F, xlab = '', ylab = '')
    legend(x = position, inset = 0, horiz = F,   # position
           col = color,                          # apply color palette
           lwd = 2,                              # apply line thicknesses
           lty = line_type,                      # apply line patterns
           title = legend_title,                 # add legend title
           legend = legend_text,                 # add legend labels
           seg.len = 2.5,                          # adjust length of lines
           cex = leg,                            # adjust legend text size
           bty = 'n')
  }
    
  
  ####
  #####
  ##### Plot age structure over time after reserve implementation #####
  N_sample <- sims_N[, , , , , indices]
  ages <- Rec_age:Max_age
  n <- length(ages)
  
  # initialize median age structures
  control_rules <- c(2, 5)
  N_medians <- array(rep(0, n*3*(Time2 + 1)*length(control_rules)), 
                     c(n, 3, Time2 + 1, length(control_rules)))
  
  # fill in median age structure
  for (i in 1:n) {
    for (a in 1:3) {
      for (t in 1:(Time2 + 1)) {
        for (cr in 1:length(control_rules)) {
          N_medians[i, a, t, cr] <- median(N_sample[i, a, Time1 + t - 1, control_rules[cr], ])
        }
      }     
    }
  }
  
  ##### Plotting parameters (age structure) #####
  
  # SAD
  SAD <- N_sample[, 1, 1, 1, 1]
  years <- seq(20, 0, by = -5)
  
  # set main titles
  CRules <- ifelse(cr == 2, 'Static', 'Transient')
  main_title1 <- paste(Species, ', ', CRules, 'Control rule \n Target DR = ', 
                       Final_DR, sep = '')
  main_title2 <- paste(Species, ', Target DR = ', Final_DR, sep = '')
  
  # use red-blue color palettes
  palette <- colorRampPalette(c('red', 'blue'))
  color1 <- palette(length(years) + 1)
  color2 <- palette(length(years))
  
  # set line types - solid for correct M, dashed for high M, dotted for low M
  line_type1 <- 1:(length(years) + 1)
  line_type2 <- 1:length(years)
  
  # set legend title and text and position
  legend_title <- expression(bold('Control Rule'))
  legend_text1  <- c("SAD", paste('Year', years))
  legend_text2  <- paste('Year', years)
  position <- 'topright'
  
  # plot margins
  plot_margins1 = c(4.5, 5, 6, 1)
  plot_margins2 = c(4.5, 6.25, 4, 1)
  
  # text sizes
  mt <- 1.75    # main title for main plot
  lab <- 1.5    # axis labels
  ax <- 1.25    # axis tick labels
  leg <- 1.25   # legend text
  
  
  ##### TODO: troubleshoot plotting one for each area for age distributions
  
  ##### Plot age distributions separately with SAD #####
  
  for (a in 1:3) {
    
    area <- ifelse(a < 2, 'far from', ifelse(a == 3, 'in', 'near'))
    sub_title <- sprintf("%s reserve", area)
    
    # y axis margins
    yy1 <- 0
    high_point <- max(N_medians[1, a, , ])
    yy2 <- round_any(high_point, 5000, f = ceiling)
    yy_by <- yy2 / 2
    
    # x axis margins
    xx1 <- Rec_age
    xx2 <- Max_age
    xx_by <- round((xx2 - xx1)/4, 0)
    
    for (cr in 1:length(control_rules)) {
      
      # set main title
      CRules <- ifelse(cr == 2, 'Static', 'Transient')
      main_title1 <- paste(Species, ', ', CRules, 
                           ' control rule \n Target DR = ', Final_DR, sep = '')
      
      # empty plot + axes + labels + titles
      par(mar = plot_margins1)
      plot(1, type = 'l',                          # make an empty line graph
           main = main_title1,                      # title of plot
           ylab = 'Median Numbers at Age', 
           xlab = 'Age',
           yaxt = 'n', xaxt = 'n',                 # get rid of y-axis
           xlim = c(xx1, xx2),                       # set x-axis limits
           ylim = c(yy1, yy2), 
           cex.main = mt, cex.lab = lab, cex.axis = ax)
      
      # set specific y-axis
      yytick <- seq(yy1, yy2, by = yy_by)              # set yaxis tick marks
      axis(side = 2,                               # specify y axis
           at = yytick,                             # apply tick marks
           labels = T,                             # apply appropriate labels
           las = 0,                                # set text horizontal
           cex.axis = ax)
      
      # set specific y-axis
      xxtick <- seq(xx1, xx2, by = xx_by)              # set yaxis tick marks
      axis(side = 1,                               # specify y axis
           at = xxtick,                             # apply tick marks
           labels = T,                             # apply appropriate labels
           las = 0,                                # set text horizontal
           cex.axis = ax) # set axis text size 
      
      # plot SAD, then add lines for 5, 10, 15, 20 years
      lines(ages, SAD, col = color1[1], lty = line_type1[1], lwd = 2)
      
      for (y in 1:length(years)) {
        lines(ages, N_medians[, a, years[y] + 1, cr], 
              lty = line_type1[y + 1], 
              col = color1[y + 1], 
              lwd = 2)
      }
      
      legend(x = position, inset = 0.01, horiz = F,   # position
             col = color1,                          # apply color palette
             lwd = 2,                              # apply line thicknesses
             lty = line_type1,                      # apply line patterns
             title = legend_title,                 # add legend title
             legend = legend_text1,                 # add legend labels
             seg.len = 3,                          # adjust length of lines
             cex = leg,                            # adjust legend text size
             bty = 'n') 
      
    }
    
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
      x2 <- Time2
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
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
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
      x2 <- Time2
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
      lines(0:Time2, rep(1, Time2 + 1), col = 'gray', lty = 3, lwd = 2)
      
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
  
  

  plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
  plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
  
  file.copy(from=plots.png.paths, 
            to="C:/Users/Vic/Google Drive/OSU/Thesis/figures/troubleshooting/no_stochasticity")
  