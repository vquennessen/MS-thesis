source('./plot_stuff.R')

filepath1 = "../../data/1e4_NM_sims_yield.Rda"
filepath2 = "../../data/1e4_NM_sims_biomass.Rda"
filepath3 = "../../data/1e4_NM_sims_SSB.Rda"

num_sims <- 1e4
A <- 5
time2 <- 20
CR <- 6

plot_stuff(filepath1, filepath2, filepath3, A, time2, CR, num_sims, 
           sample_size = num_sims, PD = 0.25)
