#Simple model of evolution of herbicide resistance as a quantitative trait. Non-spatial, annual timestep with only non-target site resistance.  
working_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_model' 
setwd(working_loc)
source('herb_resist_proccess_functions.R')
source('nonpatial_dynamic_program.R')
test_functions_broken(working_loc)
#set up evaluation points on g
eval_object = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = 0.5, seed_expantion = 5)
dg = eval_object$seed[2] - eval_object$seed[1]
## parameters
num_iter = 100 #run simulation for 50 iterations
initial_seedbank = dnorm(eval_object$seed, 0, 2) #set the intial population as a normal distribution with mean resistance value of 0 
herb_schedual = c(rep(0, 10), rep(1, 50), rep(0, 40))  #vector of 0,1's of length num_iter that defines when herbicide is applied
germination = 0.8 #germination probability
seed_survival = 0.5 #probability that a seed in the seed bank survives one year
fec_max = 100 #the maximum number of seeds per mothers
fec0 = 0 #cost of resistance when g = 0, in logits
fec_cost = 0.1 #reduction in fecundity each additional unit of g causes, in logits
additive_variance = 0.5 #variance of conditional offspring distribution
sur0 = 5 #susrvial rate when g is 0 and there is no herbicide (in logits)
sur_cost_resist = 0 #cost of higher resistenace score in terms of reduced survival (in logits)
herb_effect = 5 #effect of herbicide on survival (in logits)
survive_resist = 0.8 #protective effect of a one unit increase in resistance score g
max_sur = 0.95 #maximum survival possible
ceiling_pop = 1000000 #maximum above ground population
density_cutoff = 0.00001 #population density (on distrbution over g) above which seed evaluation points are retained in the above ground evaluation points

result = multi_iteration(num_iter = num_iter, initial_seedbank = initial_seedbank, germination = germination, eval_object = eval_object, herb_schedual = herb_schedual, sur0 = sur0, sur_cost_resist = sur_cost_resist,
  herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop, seed_survival = seed_survival, fec_max = fec_max, fec0 = fec0, 
  fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff)

seedbank_animator(results_matrix = result, eval_object = eval_object, herb_schedual = herb_schedual, pause = 0.25, xlim = c(-10, 10))




## try a simple dynamic program on IPM aimed at reducing resistance, the optimal solution will be found using Dijkstra's algorithm where the poroblem is writtena as a decision tree
## the nodes are action-timestep pairs and distance between nodes is the reward for that action timestep pair given the inital state and all the previous visited nodes 
## (i.e. actions), we want to find the biggest cumulative discounted reward. Is similar to Bellman equation but is solved foward in time rather than backwards due to difficulty in calculating 
## what the final state will be. 

##DEFINE PARAMETERS
herb_cost = 0.1
weed_impact = 0.00001
base_profit = 1 #profit from wheat with no herb and no black grass
alt_profit = 0.5
dis_rate = 0.05
time_horizon = 10
start_pop = 100
germination = 0.8 #germination probability
seed_survival = 0.5 #probability that a seed in the seed bank survives one year
fec_max = 100 #the maximum number of seeds per mothers
fec0 = 0 #cost of resistance when g = 0, in logits
fec_cost = 0.1 #reduction in fecundity each additional unit of g causes, in logits
additive_variance = 0.5 #variance of conditional offspring distribution
sur0 = 5 #susrvial rate when g is 0 and there is no herbicide (in logits)
sur_cost_resist = 0 #cost of higher resistenace score in terms of reduced survival (in logits)
herb_effect = 5 #effect of herbicide on survival (in logits)
survive_resist = 0.8 #protective effect of a one unit increase in resistance score g
max_sur = 0.95 #maximum survival possible
ceiling_pop = 1000000 #maximum above ground population
density_cutoff = 0.00001 #population density (on distrbution over g) above which seed evaluation points are retained in the above ground evaluation points

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
















