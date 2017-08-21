# plotting script for the BG population model 
# takes .csv files produced by the model, manipulates them and plots them 
# using dplyr and ggplot2

library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)


data_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output'
output_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/'

setwd(data_loc)
par_sweep = read.csv('limited_par_sweep_nonspace.csv', header = TRUE, stringsAsFactors = FALSE)

# turn the dataframe into long form
df_long = gather(par_sweep, t_lab, metric, t1:t100)

df_long = mutate(df_long, ts = as.numeric(sapply(strsplit(t_lab, 't'), FUN = function(x) x[2])),
  par_comb = paste0('fec0 = ', fec0, ' | fec_cost = ', fec_cost, ' | h_eff = ', herb_effect,
  ' | g_pro = ', g_pro, ' | os_sd = ', off_sd))

# take the data for each different measurments 
df_sur_rr = filter(df_long, measure == 'sur_rr')
df_pop = filter(df_long, measure == 'pop_size' | measure == 'ab_sur_pop')

ggplot(df_sur_rr, aes(x = ts, y = metric)) + geom_line() + 
  facet_wrap(~ par_comb)

# look at only the parameter space that is h_eff = 16
df_sur_rr_16 = filter(df_sur_rr, herb_effect == 16)

ggplot(df_sur_rr_16, aes(x = ts, y = metric)) + geom_line() + 
  facet_wrap(~ par_comb, ncol = 3)

# look at the populations 
ggplot(df_pop, aes(x = ts, y = metric, colour = measure)) + geom_line() +
 facet_wrap(~ par_comb, ncol = 3)

# look at just h_eff = 16
df_pop_16 = filter(df_pop, herb_effect == 16)

ggplot(df_pop_16, aes(x = ts, y = metric, colour = measure)) + geom_line() +
 facet_wrap(~ par_comb, ncol = 3)

# reduced plot with just g_pro = 1 | 1.5 and fec_cost = 0.35 | 0.45
df_sur = filter(df_long, measure == 'sur_rr', herb_effect == 16, fec0 == 4,
  g_pro == 1 | g_pro == 1.5, fec_cost == 0.35 | fec_cost == 0.45)
  
df_pop = filter(df_long, measure == 'pop_size' | measure == 'ab_sur_pop', herb_effect == 16, fec0 == 4,
  g_pro == 1 | g_pro == 1.5, fec_cost == 0.35 | fec_cost == 0.45)

  
sur_plt = ggplot(df_sur, aes(x = ts, y = metric)) + 
  geom_line(aes(linetype = as.factor(off_sd))) +
  facet_grid(fec_cost ~ g_pro)

pop_plt = ggplot(df_pop, aes(x = ts, y = metric, colour = measure)) + 
  geom_line(aes(linetype = as.factor(off_sd))) +
  facet_grid(fec_cost ~ g_pro)
