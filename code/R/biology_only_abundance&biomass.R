par(mfrow = c(1, 2))

# plot abundance over time 
plot(1:time, abundance_all[1, ]/1000, pch = 16, col = "blue", 
     xlab = 'Time (years)', ylab = 'Abundance (1000s of individuals)',
     yaxt = 'n', xaxt = 'n')
axis(1, seq(0, 50, 25))
axis(2, seq(0, 1250, 250))
box()

# plot biomass over time for 5 areas
plot(1:time, biomass[1, ]/1000, type = 'l', lwd = 2, col = "blue",
     xlab = 'Time (years)', ylab = 'Biomass (metric tons)', 
     yaxt = 'n', ylim = c(0, 1250), xaxt = 'n')
axis(1, seq(0, 50, 25))
axis(2, seq(0, 1250, 250))
box()