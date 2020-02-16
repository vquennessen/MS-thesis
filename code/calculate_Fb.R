devtools::install_git('https://github.com/vquennessen/densityratio.git', 
                      method = 'libcurl', force = TRUE)

library(densityratio)

# black rockfish (OR) 2015
historical_FM(Species = 'BR_OR_2015')

# cabezon (OR) 2019
historical_FM(Species = 'CAB_OR_2019')

# lingcod (OR and WA) 2017
historical_FM(Species = 'LING_OW_2017')

# canary rockfish (OR) 2015
historical_FM(Species = 'CR_OR_2015')

# english sole (Northern CA, OR, and WA) 2013
historical_FM(Species = 'ES_COW_2013')