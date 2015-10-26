#Simple model of evolution of herbicide resistance as a quantitative trait. Non-spatial, annual timestep with only non-target site resistance.  
working_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_model' 
setwd(working_loc)
source('herb_resist_proccess_functions.R')
test_functions_broken(working_loc)
#set up evaluation points on g
eval_object = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = 0.5, seed_expantion = 5)

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
ceiling_pop = 10000 #maximum above ground population
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
base_profit = 1 #profit from wheat with no herb and no black grass
alt_profit = 0.5



##DEFINE ACTIONS
# there are three actions herb0 (whet no herbicide), herb (herbicide, wheat) or alt (alternative crop which is assumed less profitable but is unaffected by black grass.)
A_herb = [0, 1, 0] # herbicide for each action
A_cost = [0, cost_herb, 0] #assume all cost and profit defined in reward function for action 3 (alt)
A_sur = [max_sur, max_sur, 0] #no above ground individuals survive to reproduction in new crop

































