SAD <- c(19849.7223, 15002.0974, 11338.3413, 8569.3340, 6476.563, 4894.8812, 
         3699.4717, 2796.0005, 2113.1717, 1597.1008, 1207.0628, 912.2785, 
         689.4852, 521.1017, 393.8402, 297.6580, 224.9651)

prop_SAD <- SAD/sum(SAD)

y1 <- 0
y2 <- 0.3
barplot(prop_SAD, main = 'Unfished Age Distribution', 
        xlab = 'Age (years)', ylab = 'Proportion at Age', 
        ylim = c(y1, y2), col = 'turquoise3', yaxt = 'n', 
        names.arg = c(as.character(1:17)), cex.names = 0.75)
axis(side = 2, at = seq(y1, y2, by = y2/2), labels = T)
box()

SAD_fished <- SAD - 1000
SAD_fished[which(SAD_fished < 0)] <- 0
prop_SAD_fished <- SAD_fished/sum(SAD_fished)

barplot(prop_SAD_fished, main = 'Fished Age Distribution', 
        xlab = 'Age (years)', ylab = 'Proportion at Age', 
        ylim = c(y1, y2), col = 'turquoise3', yaxt = 'n', 
        names.arg = c(as.character(1:17)), cex.names = 0.75)
axis(side = 2, at = seq(y1, y2, by = y2/2), labels = T)
box()
