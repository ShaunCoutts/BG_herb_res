## TODO: ACCOCIATE A POP SD WITH EACH MEAN G, FOR EACH PARAMETER SET, SINCE SD CHANGES A BIT WHEN POP STRADELS SELECTION GRADIEN 
## CAN DO THIS USING FUNCTION IN herb_resist_proccess_functions.R IN /BLACKGRASS_POPULATION_MODEL
## WOULD MAKE SENSE TO DO THIS AS PART OF THE BURNIN, AFTER THE SAMPLE POINTS ARE MADE, CALC A EQUIL SD FOR EACH SAMPL POINT
## UNDER NO_HERB AND HERB.


#Simple model of evolution of herbicide resistance as a quantitative trait. Non-spatial, annual timestep with only non-target site resistance.  
working_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation' 
setwd(working_loc)
library(colorspace)
source('herb_resist_proccess_functions_IWM.R')
source('non-spatial_dynamic_program.R')
source('IMW_DP_plotting.R')
source('non-spatial_dynamic_program_testing_diagnostics.R')
test_functions_broken(working_loc)
#set up evaluation points on g
eval_object = eval_points_builder(lower_eval_point = -20, upper_eval_point = 20, resolution = 0.5, seed_expantion = 3)
eval_object_int = eval_object
dg = eval_object$seed[2] - eval_object$seed[1]
inital_state = 100 * dnorm(eval_object$seed, 0, 0.98)
##DEFINE PARAMETERS
seed_survival = 0.3 #probability that a seed in the seed bank survives one year
fec_max = 70 #the maximum number of seeds per mothers
fec0 = 5 #cost of resistance when g = 0, in logits
fec_cost = 0.3 #reduction in fecundity each additional unit of g causes, in logits
offspring_sd = 0.7 #variance of conditional offspring distribution
sur0 = 5 #susrvial rate when g is 0 and there is no herbicide (in logits)
germination = 0.8 #germination rate
sur_cost_resist = 0 #cost of higher resistenace score in terms of reduced survival (in logits)
survive_resist = 0.5 #protective effect of a one unit increase in resistance score g
max_sur = 0.95 #maximum survival possible
pro_exposed = 0.95
dense_depend_fec = 0.0005
burnin = 100
#managment parameters
cost_plow = 10
cost_herb = 1
cost_dens = 10
mech_cost0 = 0.1
mech_cost = 0.02

income_wheat = 100
income_alt = 80
income_fallow = 0

yeild_loss_wheat = 0.01
yeild_loss_alt = 0.0001

cost_alt = 1
cost_fallow = 1

effect_herb = 10 #effect of herbicide on survival (in logits)
crop_sur_alt = 0.5
crop_sur_fal = 0
crop_fec_alt = 0.95
dens_fec_effect = 0.9
mech_sur_effect = 0.2
plow_effect = 0.3

discount_factor = 0.95
#actions space and sequences
sub_action = list(herb = c(0, 1), 
		  crop_sur = c(1, crop_sur_alt, crop_sur_fal), 
		  mech_sur = c(1, mech_sur_effect), 
		  dens_fec = c(1, dens_fec_effect), 
		  crop_fec = c(1, crop_fec_alt, 1), 
		  plow = c(0, plow_effect))
		  

#managment cost vectors
income0 = c(income_wheat, income_alt, income_fallow)
yeild_loss = c(yeild_loss_wheat, yeild_loss_alt, 0)

cost_space_nomech = list(herb = c(0, cost_herb),
			 crop = c(0, cost_alt, cost_fallow),
			 plow = c(0, cost_plow),
			 dens = c(0, cost_dens))
 
action_space = cbind(herb = c(rep(c(1, 2), each = 8), rep(c(1, 2), each = 4), 1), 
		     crop = c(rep(1, 16), rep(2, 8), 3), 
		     mech = c(rep(rep(c(1, 2), each = 2), 4), rep(c(1, 2), 4), 1), 
		     plow = c(rep(rep(c(1, 2), each = 4), 2), rep(rep(c(1, 2), each = 2), 2), 1), 
		     dens = c(rep(c(1, 2), 8), rep(1, 9)))
 
 output_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation'

# Calls the dynamic program, works much better with BRT, to run with BRT instead of gam pass the parameter BRT_para_list = list(num_trees = , int_depth = 4, learn_rate =) 
DP_policy = IMW_dynamic_program(inital_state = inital_state, germination = germination, seed_survival = seed_survival, eval_object_int = eval_object_int, pro_exposed = pro_exposed, sur0 = sur0, sur_cost_resist = sur_cost_resist, 
  effect_herb = effect_herb, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, 
  dense_depend_fec = dense_depend_fec, dg = dg, burnin = burnin, sub_action = sub_action, action_space = action_space, time_horizon = 20, discount_factor = discount_factor, 
  income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, num_samples = 500, 
  burnin_test_out_loc = output_loc, burnin_test_out_name = 'burnin_test_output.pdf', BRT_para_list = list(num_trees = 5000, int_depth = 5, learn_rate = 0.05))
 
plot_policy(Q = DP_policy$DP, output_loc = output_loc, output_name = 'policy_printout_BRT.pdf')

#make prettier plot of the policy
named_col_pal = rainbow_hcl(n = 25, start = 0, end = 300)
names(named_col_pal) <- paste0('a_', 1:25)
pretty_policy_plot(Q = DP_policy$DP[[20]], res = 100, named_col_pal = named_col_pal, output_loc = output_loc, output_name = 'policy_pretty_output.pdf', 
  xlab = 'mean g', ylab = 'population', cex.lab = 1.5)
  
#set up a plot of the optimal sub-actions across state space
act_cols = rainbow_hcl(n = 2, start = 270, end = 180)
#generate the colour pallet for the policy break down plots 
breakdown_col_pal = list(crop = rainbow_hcl(n = 3, start = 0, end = 0.66667 * 360),
  herb = act_cols, 
  plow = act_cols, 
  mech = act_cols,
  dens = act_cols) 

names(breakdown_col_pal$crop)<- c('wheat', 'alt', 'fallow') 
names(breakdown_col_pal$herb)<- c('no_herb', 'herb') 
names(breakdown_col_pal$plow)<- c('no_plow', 'plow') 
names(breakdown_col_pal$mech)<- c('no_mech', 'mech') 
names(breakdown_col_pal$dens)<- c('stand', 'high') 

#policy break down for myopic decision and long range decision
policy_breakdown_plot(Q = DP_policy$DP[[1]], action_space = action_space, res = 100, breakdown_col_pal = breakdown_col_pal, output_loc = output_loc, 
  output_name = 'policy_breakdown_output_myopic.pdf', xlab = 'mean resistance', ylab = 'population')

policy_breakdown_plot(Q = DP_policy$DP[[20]], action_space = action_space, res = 100, breakdown_col_pal = breakdown_col_pal, output_loc = output_loc, 
  output_name = 'policy_breakdown_output_T20.pdf', xlab = 'mean resistance', ylab = 'population')
  
#simulate a policy to see how the optimal policy plays out over time
sim_start_points = data.frame(mean_g = c(0, 0, 13, 13), pop = c(100, 25000, 100, 25000))
simulated_policy = policy_simulator(policy = DP_policy$DP[[20]], model_pars = DP_policy$parameters, time_period = 10, start_points = sim_start_points)
#plot the simulated policy for T = 20
plot_simulated_policy(Q = DP_policy$DP[[20]], sim_obj = simulated_policy, sims_to_plot = c(1), table_row_height = 0.05, output_loc = output_loc, 
  output_name = 'policy_trace_plot_T20_goodSP.pdf', xlim = c(0, 1), ylim = c(0, 1), cex.lab = 1.5)
plot_simulated_policy(Q = DP_policy$DP[[20]], sim_obj = simulated_policy, sims_to_plot = c(4), table_row_height = 0.05, output_loc = output_loc, 
  output_name = 'policy_trace_plot_T20_badSP.pdf', xlim = c(0, 1), ylim = c(0, 1), cex.lab = 1.5)
plot_simulated_policy(Q = DP_policy$DP[[20]], sim_obj = simulated_policy, sims_to_plot = c(3), table_row_height = 0.05, output_loc = output_loc, 
  output_name = 'policy_trace_plot_T20_low_pop_hr_SP.pdf', xlim = c(0, 1), ylim = c(0, 1), cex.lab = 1.5)
plot_simulated_policy(Q = DP_policy$DP[[20]], sim_obj = simulated_policy, sims_to_plot = c(2), table_row_height = 0.05, output_loc = output_loc, 
  output_name = 'policy_trace_plot_T20_high_pop_lowr_SP.pdf', xlim = c(0, 1), ylim = c(0, 1), cex.lab = 1.5)
#plot the simulated policy for myopic decision
simulated_policy = policy_simulator(policy = DP_policy$DP[[1]], model_pars = DP_policy$parameters, time_period = 10, start_points = sim_start_points)
#plot the simulated policy for T = 20
plot_simulated_policy(Q = DP_policy$DP[[1]], sim_obj = simulated_policy, sims_to_plot = c(1), table_row_height = 0.05, output_loc = output_loc, 
  output_name = 'policy_trace_plot_T1_goodSP.pdf', xlim = c(0, 1), ylim = c(0, 1), cex.lab = 1.5)
plot_simulated_policy(Q = DP_policy$DP[[1]], sim_obj = simulated_policy, sims_to_plot = c(4), table_row_height = 0.05, output_loc = output_loc, 
  output_name = 'policy_trace_plot_T1_badSP.pdf', xlim = c(0, 1), ylim = c(0, 1), cex.lab = 1.5)
plot_simulated_policy(Q = DP_policy$DP[[1]], sim_obj = simulated_policy, sims_to_plot = c(3), table_row_height = 0.05, output_loc = output_loc, 
  output_name = 'policy_trace_plot_T1_low_pop_hr_SP.pdf', xlim = c(0, 1), ylim = c(0, 1), cex.lab = 1.5)
plot_simulated_policy(Q = DP_policy$DP[[1]], sim_obj = simulated_policy, sims_to_plot = c(2), table_row_height = 0.05, output_loc = output_loc, 
  output_name = 'policy_trace_plot_T1_high_pop_lowr_SP.pdf', xlim = c(0, 1), ylim = c(0, 1), cex.lab = 1.5)

#make simple plot of value surface 
plot_value_surface(Q = DP_policy$DP[[20]], output_loc = output_loc, output_name = 'value_surface_T20.pdf', xlab = 'resistance', ylab = 'population', 
  cex.lab = 1.5, cex.main = 3, main = 'Value surface 20 year time horizon')

plot_value_surface(Q = DP_policy$DP[[1]], output_loc = output_loc, output_name = 'value_surface_T1.pdf', xlab = 'resistance', ylab = 'population', 
  cex.lab = 1.5, cex.main = 3, main = 'Value surface myopic')
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#trying to fugure out the simulated popuation going the wrong way, getting less resistant with more herbicide
#multi_iteration() produces the expected result, but the simulation DP function does the opposite, also even when there 
#is no fec cost then the line just goes stright up, not getting more resistant with each herb use. So there is no 
#force pushing the resistance up it seems. need to figure it out tomorrow.

fec_test = fecundity_closure(eval_points_object = eval_object_int, offspring_sd = offspring_sd, fec0 = fec0, fec_cost = fec_cost)
current_state = build_state_1level(dist_top = 100 * dnorm(eval_object_int$seed, -4, 0.97), eval_points = eval_object_int$seed, dg = dg)
test_action = c(herb = 2, crop = 1, mech = 1, dens = 1, plow = 1)

state_reward_next_econ_plot(current_state = current_state, eval_object = eval_object_int, action = test_action, sub_action = sub_action, seed_survival = seed_survival, 
  germination = germination, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, 
  survive_resist = survive_resist, fec_max = fec_max, dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, 
  mech_cost0 = mech_cost0, mech_cost = mech_cost, fec_function = fec_test, dg = dg)
  
#problem does not seem to be in the state update, which behaves as expected, so next possibility is the policy_simulator
test = policy_simulator_test(policy = DP_policy$DP[[10]], model_pars = DP_policy$parameters, time_period = 3, start_points = data.frame(mean_g = 0, pop = 1000), 
  output_loc = output_loc, output_name = 'policy_simulation_test.pdf')

int_pop = 100
seedbank_int = int_pop * dnorm(eval_object_int$seed, 0, 0.98)
num_iter = 50
action_seq = cbind(herb = rep(1, num_iter), crop = rep(1, num_iter), mech = rep(1, num_iter), plow = rep(1, num_iter), dens = rep(1, num_iter))
#test the intial iterations to see how the burnin pop should react under repeated herbicide or non-herbicide 
test_burnin_1 = multi_iteration_plot(seedbank_initial = seedbank_int, germination = germination, seed_survival = seed_survival, eval_object = eval_object_int, 
  pro_exposed = pro_exposed, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, max_sur = max_sur, 
  fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, dg = dg, num_iter = num_iter, 
  sub_action = sub_action, action_seq = action_seq, output_loc = output_loc, output_name = 'multi_iteration_test_output.pdf')
  
action_seq = cbind(herb = rep(2, num_iter), crop = rep(1, num_iter), mech = rep(1, num_iter), plow = rep(1, num_iter), dens = rep(1, num_iter))
test_burnin = multi_iteration(seedbank_initial = seedbank_int, germination = germination, seed_survival = seed_survival, eval_object = eval_object_int, 
  pro_exposed = pro_exposed, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, max_sur = max_sur, 
  fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, dg = dg, num_iter = num_iter, 
  sub_action = sub_action, action_seq = action_seq)
  
rowSums(test_burnin_1)*dg
rowSums(test_burnin)*dg

#NOTE: looks like the way density dependence is modelled is very important since can get big over-yeilding situations where applying herbicide 
# makes the seedbank much worse? is this realistic? depends very strongly on seeds per plant 

#NOTES ON COSTINg
#Fixed cost: machine cost -> depreciated and intial cost, pay no matter what, this might fall out as relative cost that matters
#varible costs: labour, chemicals, fule,look for tractor hours -> cerials -> 9hr per ha -> translate that to money
#machines: tractor and combine harvester

#use gross margins of alt crop vs wheat, based on average yeild, average price, nix has gross margin. Also cropping costs like fungiciced
#also fallow comes with subisdies, low land farmrs £207/ha from nix
#red desile £0.8 per liter

#discount rate, can include inflation and interest rates, UK inflation 2.5%