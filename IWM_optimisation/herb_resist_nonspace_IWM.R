#Simple model of evolution of herbicide resistance as a quantitative trait. Non-spatial, annual timestep with only non-target site resistance.  
working_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation' 
setwd(working_loc)
source('herb_resist_proccess_functions_IWM.R')
source('non-spatial_dynamic_program.R')
source('non-spatial_dynamic_program_testing_diagnostics.R')
test_functions_broken(working_loc)
#set up evaluation points on g
eval_object = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = 0.5, seed_expantion = 5)
eval_object_int = eval_object
dg = eval_object$seed[2] - eval_object$seed[1]
inital_state = 100 * dnorm(eval_object$seed, 0, 2)
##DEFINE PARAMETERS
weed_impact = 0.00001
dis_rate = 0.05
time_horizon = 10
start_pop = 100
seed_survival = 0.5 #probability that a seed in the seed bank survives one year
fec_max = 100 #the maximum number of seeds per mothers
fec0 = 0 #cost of resistance when g = 0, in logits
fec_cost = 0.1 #reduction in fecundity each additional unit of g causes, in logits
offspring_sd = 0.7 #variance of conditional offspring distribution
sur0 = 5 #susrvial rate when g is 0 and there is no herbicide (in logits)
germination = 0.8 #germination rate
sur_cost_resist = 0 #cost of higher resistenace score in terms of reduced survival (in logits)
survive_resist = 0.8 #protective effect of a one unit increase in resistance score g
max_sur = 0.95 #maximum survival possible
density_cutoff = 0.00001 #population density (on distrbution over g) above which seed evaluation points are retained in the above ground evaluation points
seed_movement = 0.8
pro_exposed = 0.8
dense_depend_fec = 0.002
burnin = 20
#managment parameters
effect_plow = 0.9 #proprtion of seed moved from one depth level to another by plowing 
effect_herb = 5 #effect of herbicide on survival (in logits)

cost_plow = 0.05
cost_herb = 0.1
cost_dens = 0.5
mech_cost0 = 0.01
mech_cost = 0.05

income_wheat = 100
income_alt = 50
income_fallow = 0

yeild_loss_wheat = 0.0001
yeild_loss_alt = 0.000001

cost_alt = 0.5
cost_fallow = 0.5

crop_sur_alt = 0.8
crop_sur_fal = 0
crop_fec_alt = 0.95
dens_fec_effect = 0.9
mech_sur_effect = 0.1
plow_effect = 0.3

discount_factor = 0.95
#actions space and sequences
sub_action = list(herb = c(0, 1), 
		  crop_sur = c(1, crop_sur_alt, crop_sur_fal), 
		  mech_sur = c(1, mech_sur_effect), 
		  dens_fec = c(1, dens_fec_effect), 
		  crop_fec = c(1, crop_fec_alt, 1), 
		  plow = c(0, plow_effect))
		  

action_seq = cbind(herb = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2), 
		   crop = c(1, 1, 1, 2, 2, 3, 1, 1, 1, 1), 
		   mech = c(1, 2, 1, 1, 1, 1, 1, 1, 1, 2), 
		   dens = c(1, 1, 2, 1, 1, 1, 2, 2, 2, 1), 
		   plow = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2))	  
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

# plot of the seed bank through time with lots of differnt plots to see what is happening behind the scenes 
single_iteration_plot(seedbank_current = seedbank_current, germination = germ_rate_wheat, mech_control = 1, crop_effect_sur = 1, seed_survival = seed_survival, 
  seed_movement = seed_movement, eval_object = eval_object, pro_exposed = pro_exposed, herb_rate = 1, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, 
  survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, 
  density_cutoff = density_cutoff, crop_effect_fec = 1, density_effect_fec = 1, dg = dg)
 
multi_iteration_plot(seedbank_initial = seedbank_current, germination = germination, seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, sur0 = sur0, 
  sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, 
  offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, density_cutoff = density_cutoff, dg = dg, num_iter = dim(action_seq)[1], sub_action = sub_action, 
  action_seq = action_seq, output_loc = output_loc, output_name = 'multi_iteration_output.pdf')
 
IMW_dynamic_program(inital_state = inital_state, germination = germination, seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, sur0 = sur0, 
  sur_cost_resist = sur_cost_resist, effect_herb = effect_herb, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, 
  offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, density_cutoff = density_cutoff, dg = dg, burnin = burnin, sub_action = sub_action, 
  noherb_action_seq = noherb_action_seq, burnin_test_out_loc = output_loc, burnin_test_out_name = 'burnin_test_output.pdf')
 








##DEFINE PARAMETERS
herb_cost = 0.1
weed_impact = 0.00001
dis_rate = 0.05
time_horizon = 10
start_pop = 100
seed_survival = 0.5 #probability that a seed in the seed bank survives one year
fec_max = 100 #the maximum number of seeds per mothers
fec0 = 0 #cost of resistance when g = 0, in logits
fec_cost = 0.1 #reduction in fecundity each additional unit of g causes, in logits
additive_variance = 0.5 #variance of conditional offspring distribution
sur0 = 5 #susrvial rate when g is 0 and there is no herbicide (in logits)
sur_cost_resist = 0 #cost of higher resistenace score in terms of reduced survival (in logits)
survive_resist = 0.8 #protective effect of a one unit increase in resistance score g
max_sur = 0.95 #maximum survival possible
ceiling_pop = 1000000 #maximum above ground population
density_cutoff = 0.00001 #population density (on distrbution over g) above which seed evaluation points are retained in the above ground evaluation points

#managment parameters
germ_rate_wheat = 0.8
germ_rate_alt = 0.1
germ_rate_fallow = 0

effect_plow = 0.9 #proprtion of seed moved from one depth level to another by plowing 
effect_herb = 5 #effect of herbicide on survival (in logits)

cost_plow = 0.05
cost_herb = 0.1

income_wheat = 10
income_alt = 2
income_fallow = 0



## define the action space a three by 10 matrix with 3 sub-actions all combinations except 2
#sub action effects
sub_action_germ = c(germ_rate_wheat, germ_rate_alt, germ_rate_fallow)
sub_action_plow = c(0, effect_plow)
sub_action_herb = c(0, 1)
#A is the action space and calls the sub-actions by index to the sub actions vectors the order is crop|herb|plow 
A = rbind(c(1, 1, 1),
	  c(1, 1, 2), 
	  c(1, 2, 1),
	  c(1, 2, 2),
	  c(2, 1, 1),
	  c(2, 1, 2),
	  c(2, 2, 1),
	  c(2, 2, 2),
	  c(3, 2, 2))



## try a simple dynamic program on IPM aimed at reducing resistance, the optimal solution will be found using Dijkstra's algorithm where the poroblem is writtena as a decision tree
## the nodes are action-timestep pairs and distance between nodes is the reward for that action timestep pair given the inital state and all the previous visited nodes 
## (i.e. actions), we want to find the biggest cumulative discounted reward. Is similar to Bellman equation but is solved foward in time rather than backwards due to difficulty in calculating 
## what the final state will be. 


##DEFINE ACTIONS
# there are three actions herb0 (whet no herbicide), herb (herbicide, wheat) or alt (alternative crop which is assumed less profitable but is unaffected by black grass.)
A_herb = c(0, 1, 0) # herbicide for each action
A_cost = c(0, herb_cost, 0) #assume all cost and profit defined in reward function for action 3 (alt)
A_sur = c(max_sur, max_sur, 0) #no above ground individuals survive to reproduction in new crop

#start by letting the population settle to an equilibrium state 
inital_pop = multi_iteration(num_iter = 20, initial_seedbank = dnorm(eval_object$seed, 0, 2), germination = germination, eval_object = eval_object, herb_schedual = rep(0, 20), 
  sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop, seed_survival = seed_survival, 
  fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff)
inital_pop = inital_pop[dim(inital_pop)[1], ]
#reduce pop so it is at a set size
inital_pop = inital_pop * (start_pop / sum(inital_pop))

#start building the tree
node_list = list()
node_list[[1]] = intial_node(state0 = inital_pop)

#make a list to hold any nodes that get to the time_horizon
reached_horizon = list()
min_val_at_hor = Inf

cont_search = TRUE
next_source = 1
#search the tree for the best path to the low value at time_horizon
while(cont_search){
  print(length(node_list))
  #add three child nodes, one for each action
  for(i in seq_along(A_cost)){
    node_list <- add_node(node_list = node_list, act = i, parent_index = next_source, alt_profit = alt_profit, seed_survival = seed_survival, base_profit = base_profit, 
      weed_impact = weed_impact, germination = germination, eval_object = eval_object, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, 
      max_sur = max_sur, ceiling_pop = ceiling_pop, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff, 
      dg = dg, dis_rate = dis_rate)
  }
  node_list = node_list[-next_source] #remove the visited node from the list
  new_values_set = sapply(node_list, FUN = function(x) x$value)
  next_source = which(new_values_set == min(new_values_set))
  #check if we have reached time horizon and if have find the next source that has not reached the time_horizon
  while(node_list[[next_source]]$time >= time_horizon){
    if(node_list[[next_source]]$value < min_val_at_hor) min_val_at_hor = node_list[[next_source]]$value
    reached_horizon[[length(reached_horizon) + 1]] = node_list[[next_source]]
    node_list = node_list[-next_source]
    if(length(node_list) > 0){
      new_values_set = sapply(node_list, FUN = function(x) x$value)
      next_source = which(new_values_set == min(new_values_set))
    }else{
      cont_search = FALSE
      break
    }
  }
  if(length(node_list) == 0) break
  #check if the next source has a value lower than the min value at time_horizon (we know these branches are dead ends), find next viable source if it exists
  while(node_list[[next_source]]$value > min_val_at_hor){
    node_list = node_list[-next_source] #remove the dead end from the search
    if(length(node_list) > 0){
      new_values_set = sapply(node_list, FUN = function(x) x$value)
      next_source = which(new_values_set == min(new_values_set))
    }else{
      cont_search = FALSE
      break
    }
  }
}

#find the best one at the time horizion (will normally only be one the reaches that far)
best_solution = reached_horizon[[which(sapply(reached_horizon, FUN = function(x) x$value) == min(sapply(reached_horizon, FUN = function(x) x$value)))]]
action_sequence = rev(best_solution$acts)
#simulate the best solution
solution_nodes = list()
solution_nodes[[1]] = intial_node(state0 = inital_pop) 
for(i in 2:length(action_sequence)){
  solution_nodes = add_node(node_list = solution_nodes, act = action_sequence[i], parent_index = i - 1, alt_profit = alt_profit, seed_survival = seed_survival, base_profit = base_profit, weed_impact = weed_impact, 
      germination = germination, eval_object = eval_object, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, 
      max_sur = max_sur, ceiling_pop = ceiling_pop, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff, 
      dg = dg, dis_rate = dis_rate) 
}
#animated visulisation
simulate_solution(solution_nodes, eval_object, pause = 0.5)
















