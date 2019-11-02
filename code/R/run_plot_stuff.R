source('./plot_stuff.R')

filepath1 = "../../data/1e2_NM_sims_yield.Rda"
filepath2 = "../../data/1e2_NM_sims_biomass.Rda"
filepath3 = "../../data/1e2_NM_sims_SSB.Rda"
num_sims <- 1e2

plot_stuff(filepath1, filepath2, filepath3, A, time2, CR, num_sims, 
           sample_size = num_sims, PD = 0.25)