# plotting script for the BG population model 
# takes .csv files produced by the model, manipulates them and plots them 
# using dplyr and ggplot2

library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

data_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output'
output_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/'

# helper function to get survival from g
g_2_sur = function(g, rho, s0, hef){

  return(1 / (1 + exp(-(s0 - (hef - g * rho)))))

}


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


  
# plot of survival, seed bank size and above ground after control population size
df_sur_pop = filter(df_long, measure == 'pop_size' | measure == 'ab_sur_pop' | measure == 'sur_rr', herb_effect == 16, 
  fec0 == 4, g_pro == 1 | g_pro == 1.5, fec_cost == 0.45)

df_sur_pop = mutate(df_sur_pop, pretty_measure = sapply(measure, function(x){
  
  if(x == 'sur_rr'){
  
    return('survival rr')
    
  }else{
    
    return('population size')
    
  }
}))
 
pop_plt = ggplot(df_sur_pop, aes(x = ts, y = metric, colour = measure)) + 
  geom_line(aes(linetype = as.factor(off_sd))) +
  facet_grid(pretty_measure ~ g_pro, scales = 'free_y')

  
## more focused run #######################################################################

setwd(data_loc)
par_sweep = read.csv('rho_Va_par_sweep_nonspace.csv', header = TRUE, stringsAsFactors = FALSE)

# turn the dataframe into long form
df_long = gather(par_sweep, t_lab, metric, t1:t100)

df_long = mutate(df_long, ts = as.numeric(sapply(strsplit(t_lab, 't'), FUN = function(x) x[2])))

df_sur_pop = filter(df_long, measure == 'pop_size' | measure == 'ab_sur_pop' | measure == 'sur_rr', 
  g_pro == 1 | g_pro == 1.5)
  
# add some pretty labels and lest the faceting work porperly
df_sur_pop = mutate(df_sur_pop, pretty_measure = as.character(sapply(measure, function(x){
  
  if(x == 'sur_rr'){
  
    return('survival rr')
    
  }else{
    
    return('population size')
    
  }
})), Va_pretty = paste0('Offspring var = ', off_sd^2),
rho_pretty = paste0('rho = ', g_pro))
 
 
pop_plt = ggplot(df_sur_pop, aes(x = ts, y = metric, colour = measure)) + 
  geom_line(aes(linetype = Va_pretty)) +
  facet_grid(pretty_measure ~ rho_pretty, scales = 'free_y')

####################################################################################################################
# look at the seed translocation experements
setwd(data_loc)
trans_expr = read.csv('rho_Va_trans_expr.csv', header = TRUE, stringsAsFactors = FALSE)

# turn the dataframe into long form
df_expr = gather(trans_expr, t_lab, metric, t1:t120)

df_expr = mutate(df_expr, ts = as.numeric(sapply(strsplit(t_lab, 't'), FUN = function(x) x[2])),
  ts_inj = ts - est_period, inj_sur_rr = g_2_sur(inj_g, g_pro, s0, herb_effect), 
  inj_rr_pretty = paste0('survival rr inj. = ', round(inj_sur_rr, 2)), inj_R_pretty = paste0('%R inj. = ', injRR / 10), 
  targ_scen = ifelse(intrr == 0, 'empty', ifelse(intrr > 0 & herb1 == 1, 'naive', 'exposed')),
  Va_pretty = paste0('Offspring var = ', off_sd ^ 2))

# select some of the parameter space to plot
df_plot = filter(df_expr, g_pro == 1.5, measure == "sur_rr" | measure == "pro_R",
  inj_R_pretty == '%R inj. = 0.1', ts >= est_period)

df_plot = mutate(df_plot, pretty_metric = ifelse(measure == 'sur_rr', 'survival rr', '%R'))
  
  
trans_expr_plt = ggplot(df_plot, aes(x = ts_inj, y = metric * 100, colour = Va_pretty)) +
  geom_line(aes(linetype = pretty_metric), size = 1.3) + labs(x = 'time step', y = 'percent (%)') + 
  annotate('text', label = paste0(letters[1:9], ')'), size = 5, x = 0, y = 100) + 
  theme(legend.position = c(0.8, 0.5), panel.background = element_rect(fill = grey(0.95)),
    panel.grid.major = element_line(colour = grey(1)),
    legend.background = element_rect(fill = grey(1))) +
  facet_grid(targ_scen ~ inj_rr_pretty)

setwd(output_loc)
pdf('TSR_rho1.5_injRR01.pdf', width = 12, height = 12)

  trans_expr_plt

dev.off()

## look at the effect of rho on evoloution of TSR in a population
# look at the seed translocation experements
setwd(data_loc)
trans_expr = read.csv('rho_Va_trans_expr.csv', header = TRUE, stringsAsFactors = FALSE)

# turn the dataframe into long form
df_expr = gather(trans_expr, t_lab, metric, t1:t120)

# take a time slice and two metrics to plot
df_rho = filter(df_expr, t_lab == paste0('t', est_period + 100),
  measure == "sur_rr" | measure == "pro_R")

df_rho = mutate(df_rho, ts = as.numeric(sapply(strsplit(t_lab, 't'), FUN = function(x) x[2])),
  ts_inj = ts - est_period, inj_sur_rr = g_2_sur(inj_g, g_pro, s0, herb_effect), 
  inj_rr_pretty = paste0('survival rr inj. = ', round(inj_sur_rr, 2)), inj_R_pretty = paste0('%R inj. = ', injRR / 10), 
  targ_scen = ifelse(intrr == 0, 'empty', ifelse(intrr > 0 & herb1 == 1, 'naive', 'exposed')),
  Va_pretty = paste0('Offspring var = ', off_sd ^ 2), pretty_metric = ifelse(measure == 'sur_rr', 'survival rr', '%R'))

df_rho = filter(df_rho, inj_R_pretty == '%R inj. = 0.1')
  
rho_plt = ggplot(df_rho, aes(x = g_pro, y = metric * 100, colour = Va_pretty)) + 
  geom_line(aes(linetype = pretty_metric), size = 1.3) +  labs(x = expression(rho), y = 'percent (%)') + 
  annotate('text', label = paste0(letters[1:9], ')'), size = 5, x = 1, y = 100) + 
  theme(legend.position = c(0.15, 0.85), panel.background = element_rect(fill = grey(0.95)),
    panel.grid.major = element_line(colour = grey(1)),
    legend.background = element_rect(fill = grey(1))) +
  facet_grid(targ_scen ~ inj_rr_pretty)

setwd(output_loc)
pdf('TSR_quant_rho_injRR01.pdf', width = 12, height = 12)

  rho_plt

dev.off()
  
## a run with no-TSR, just shows the evolution of quantitative resistance 
setwd(data_loc)
noTSR = read.csv("rho_Va_noTSR.csv", header = TRUE, stringsAsFactors = FALSE)

# turn the dataframe into long form
df_noTSR = gather(noTSR, t_lab, metric, t1:t120)

df_noTSR = mutate(df_noTSR, ts = as.numeric(sapply(strsplit(t_lab, 't'), FUN = function(x) x[2])),
  ts_inj = ts - est_period, inj_sur_rr = g_2_sur(inj_g, g_pro, s0, herb_effect), 
  inj_rr_pretty = paste0('survival rr inj. = ', round(inj_sur_rr, 2)),
  Va_pretty = paste0('Offspring var = ', off_sd ^ 2))

df_noTSR = filter(df_noTSR, g_pro == 1.5, measure == "sur_rr" | measure == "pop_size" | measure == 'ab_sur_pop',
  ts >= est_period, inj_rr_pretty == 'survival rr inj. = 0.01')

df_noTSR = mutate(df_noTSR, pretty_measure = sapply(measure, function(x){
  
  if(x == 'sur_rr'){
  
    return('survival rr')
    
  }else{
    
    return('population size')
    
  }
}),
pretty_metric = sapply(measure, function(x){
  
  if(x == 'sur_rr'){
  
    return('survival rr')
    
  }else{if(x == 'pop_size'){
  
    return('seed bank')
  
  }else{
  
    return('above ground pop.')
  
  }}
}))
 
noTSR_plt = ggplot(df_noTSR, aes(x = ts_inj, y = metric, colour = Va_pretty)) +
  geom_line(aes(linetype = pretty_metric)) +
  labs(x = 'time step', y = 'proportion surviving                                                           population') + 
  theme(legend.position = c(0.8, 0.15), panel.background = element_rect(fill = grey(0.95)),
    panel.grid.major = element_line(colour = grey(1)),
    legend.background = element_rect(fill = grey(1))) + 
  facet_grid(pretty_measure ~ ., scales = 'free_y') 
  
setwd(output_loc)
pdf('noTSR_rho15_injrrg01.pdf', width = 7, height = 10)
  noTSR_plt
dev.off()

#######################################################################################################################################################################
# spatial model
setwd(data_loc)
TSR_space = read.csv("rho_Va_TSR_space.csv", header = TRUE, stringsAsFactors = FALSE)

# find the first index that have a non-zero value of num_RR at center cell to get introduction time
inj_t = which(TSR_space$x100 > 0)[1]

# turn the dataframe into long form
df_TSR_space = gather(TSR_space, loc_lab, value, x1:x200)

# add a few pretty labels and change some coloumns to numeric
df_TSR_space = mutate(df_TSR_space, loc = as.numeric(sapply(strsplit(loc_lab, 'x'), FUN = function(x) x[2])),
  ts_inj = ts - inj_t, inj_sur_rr = g_2_sur(inj_g, g_pro, s0, herb_effect),
  inj_rr_pretty = paste0('survival rr inj. = ', round(inj_sur_rr, 2)),
  Va_pretty = paste0('Offspring var = ', off_sd ^ 2))

# define some parameter space to plot over and select some time periods to look at 
# since we are plotting over space and looking at time slices, also take a couple of
# metrics %R and sur_rr
plot_df = filter(df_TSR_space, inj_rr_pretty == 'survival rr inj. = 0.5', 
  ts_inj == 5 | ts_inj == 25 | ts_inj == 50 | ts_inj == 75 | ts_inj == 100,
  metric == 'pro_R', g_pro == 1.5)

# make the colour pallet for plot with custom limits
my_blues = brewer.pal(n = 7, "Blues")[3:7]
  
  
TSR_space_plt = ggplot(plot_df, aes(x = loc, y = value * 100, colour = as.factor(ts_inj))) + 
  geom_line() +   labs(x = 'location', y = '%R') + 
  theme(legend.position = c(0.1, 0.70), panel.background = element_rect(fill = grey(0.95)),
    panel.grid.major = element_line(colour = grey(1)),
    legend.background = element_rect(fill = grey(1))) + 
    annotate('text', label = paste0(letters[1:9], ')'), size = 5, x = 1, y = 100) + 
 scale_colour_manual(values = my_blues) + facet_grid(scen ~ Va_pretty)




























