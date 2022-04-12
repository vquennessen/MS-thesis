Time1 = 50
TimeT = 70
Final_DR = 0.6
Nat_mortality = c(0.09, 0.14, 0.19)
nm = 1

years <- 0:(TimeT - Time1)
target_DR <- 1 - (1 - Final_DR)*(1 - exp(-1 * Nat_mortality[nm] * years))

par(mar = c(4.7, 4.5, 3, 0.2))
plot(years, target_DR, main = 'Target Density Ratios',
     ylab = 'Target Density Ratio',
     xlab = 'Years since Reserve Implementation',
     pch = 15, col = 'blue', ylim = c(0.6, 1), 
     cex.lab = 1.5)
abline(h = 0.6, col = 'red', lwd = 2)
legend(x = 'topright', inset = 0.05, lwd = 2, col = c('blue', 'red'),
       c('Transient CR', 'Static CR'))
