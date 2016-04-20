# Run the population model and to see how the population behaves and if the populaiton 
# can over yeild. We also use some field data to filter out parameter combinations 
# that produce unrealistic outcomes. We then use a sensitivtiy anaylsis to determine 
# which parameters are important in driving population growth, time to K and amount of 
# over-yeilding.

#file locations
code_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model'
out_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output'

#libraries and source codes
library(lhs)
setwd(code_loc)
library(parallel)
source('herb_resist_proccess_functions.R')
source('pop_model_plot.R')


# population model parameter est and ranges
seed_sur = 0.45
germ_prob = 0.52
fec_max = 45
fec_dd = 0.004
fec0 = 8
offspring_sd = 1
sur0 = 10
pro_exposed = 0.8


param_ranges = rbind(seed_sur = c(0.22, 0.79), germ_prob = c(0.45, 0.6), fec_max = c(30, 300), 
  fec_dd = c(0.0001, 0.01), fec0 = c(5, 10), pro_exposed = c(0.5, 1), fec_cost_mult = c(0.1, 2),
  herb_effect_mult = c(2, 3), sur_protect_mult = c(0.1, 2), int_Rr = c(1.0e-06, 1.0e-01))

# model run setup 
intial_pop_size = 100
dg = 0.5
eval_object = eval_points_builder(lower_eval_point = -20, upper_eval_point = 20, resolution = dg, 
  seed_expantion = 3)
seedbank_initial = dnorm(eval_object$seed, 0, 1.391607) * intial_pop_size
int_Rr = 1.0e-06
TSR_inital = c(RR = 0, Rr = int_Rr, rr = 1 - int_Rr)
res_genotypes = c('RR')
num_iter = 100
num_par_comb = 100

pop_plotter(pop_list = pop_list, eval_points = eval_object$above_ground, dg = dg, out_loc = out_loc)
  
  
get_equil_sd(pop_size = 2000, mean_g = 0, int_sd = 1.239828, germination = germ_prob, 
  seed_survival = seed_sur, eval_object = eval_object, pro_exposed = 0.8, 
  sur0 = sur0, herb_effect = sur0 * 2.5, survive_resist = sur0 * 2.5 * 1, fec_max = fec_max, 
  fec0 = fec0, fec_cost = fec0 * 0.1, offspring_sd = 1, dense_depend_fec = fec_dd, 
  dg = dg, num_iter = num_iter)$sd_g
  
# It appears that the amount of selection the populatin is under affects the equlibrium SD   
# run each parameter combination 100 ts so the varaince can stabilise, then kill all but intial_pop_size 
# individuals and start measuring from then. check how sd reacts, may be also a funciton of populaion size 
# in that seg. var. only 2 add. var at larger pops.

#generate a set of parameters evenly sampled over the parameter space
system.time({
sampled_points = randomLHS(n = num_par_comb, k = 10)#, dup = 2)
sampled_param = lapply(1:dim(sampled_points)[1], FUN = function(x){
  ((param_ranges[, 2] - param_ranges[, 1]) * sampled_points[x, ]) + param_ranges[, 1]
}) 
})

test = mclapply(sampled_param, FUN = function(x){ 
  pop_run_TSR(param_vect = x, int_pop_size = intial_pop_size, int_sd_g = 1.4, 
  res_genotypes = res_genotypes, sur0 = 10, offspring_sd = 1, burnin_time = 100, num_iter = 50, 
  eval_object = eval_object, dg = dg)
}, mc.cores = 2, mc.allow.recursive = FALSE)

test= pop_run(param_vect = sampled_param[[1]], int_pop_size = 100, int_sd_g = 1.4, sur0 = 10, offspring_sd = 1, 
  burnin_time = 100, num_iter = 50, eval_object = eval_object, dg = dg)


# get the output object form ICEBERG
setwd(code_loc)
load('parameter_space_test_ZHUMAC.Rdata')

#find the population of each parameter set over time 
population_out = pops_list_cleanup(parameter_space_test)
pop_meseaurs = pop_metrics_simple(population_out)

# now filter out those parameter combinations that result in unrealistic populations
# using the field data to define what is unrealistic we get 
# a min of 16,348 and a max of 132,003 plants in resistant populations 
# This could be taken as the max pop under herb or the final above ground population 
# under herb  
LOW_POP = 16348
HIGH_POP = 132003
filtered_data = as.data.frame(pop_meseaurs[(pop_meseaurs[, 'fin_ag_post_herb'] > LOW_POP & pop_meseaurs[, 'fin_ag_post_herb'] < HIGH_POP) | 
  (pop_meseaurs[, 'max_ag_post_herb'] > LOW_POP & pop_meseaurs[, 'max_ag_post_herb'] < HIGH_POP), ])

#from 3000 parameter combinations only 673 managed to produce results in the very large range observed 
summary(filtered_data)

# most of these above ground population estimates come in at the low end suggesting soem of the parmeter ranges need to be
# expanded a bit to make larger populations easier to get in particular all selected populations had 
# fec_dd < 0.000875, so re-do the parameter sweep with fec_dd = [0.000001, 0.001], also maybe open up the 
# effect of herbicide to be even stronger 

filtered_inds = which((pop_meseaurs[, 'fin_ag_post_herb'] > LOW_POP & pop_meseaurs[, 'fin_ag_post_herb'] < HIGH_POP) | 
  (pop_meseaurs[, 'max_ag_post_herb'] > LOW_POP & pop_meseaurs[, 'max_ag_post_herb'] < HIGH_POP))

  
TSR_pop_v_time(pops_list = population_out, plot_inds = filtered_inds, out_loc = out_loc, out_name = 'pop_v_time_TSR_all.pdf')

 
  
  
  
  
  
  
  
  
  
  