#functions for the dynamic programing of the non-spatial model of herb-resistance 

## creat a node with all the information needed.
node <- function(act, parent, alt_profit, min_reward, seed_survival, base_profit, herb_cost, dg, dis_rate, ...){
  if(act == 3){
    new_state = seedbank_update(parent$state, seed_survival)
    value = parent$value + ((base_profit - alt_profit) / ((1 + dis_rate) ^ (parent$time_step + 1)))
  }else{
    N1 = state_transition(seedbank_current = parent$state, herb_rate = A_herb[act], seed_survival = seed_survival, ...){
    new_state N1$new_state
    value = parent$value + ((base_profit - reward(act, N1$survivors)) / ((1 + dis_rate) ^ (parent$time_step + 1)))
  }
  
  list(acts = c(act, partent$acts), time_step = parent$time_step + 1, value = value, state = new_state)
}

intial_node <- function(state0){
  list(acts = 0, time_step = 0, value = 0, state = state0)
}

reward <- function(act, survivors, base_profit, herb_cost, dg){
  N = sum(survivors) * dg
  max(0, base_profit - N * herb_cost) - A_cost[act] 
}

seedbank_update(seedbank_current, seed_survival){
  seedbank_current * seed_survival
}

state_transition <- function(seedbank_current, germination, eval_object, herb_rate, sur0, sur_cost_resist, herb_effect, survive_resist, max_sur, ceiling_pop, seed_survival, fec_max, fec0, fec_cost, additive_variance, density_cutoff){
  new_plants = emergence(seedbank_current = seedbank_current, germination = germination) 
  eval_object = eval_points_update(eval_points_object = eval_object, new_plants, density_cutoff = density_cutoff) #update evaluation window
  survivors = survival(N_0 = new_plants[eval_object$above_ground_index], eval_points = eval_object$above_ground, herb_rate = herb_rate, sur0 = sur0, sur_cost_resist = sur_cost_resist, 
	  herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop) 
  new_seedbank = seedbank(seedbank0 = seedbank_current, seed_survival = seed_survival, germination = germination, eval_object = eval_object, 
	  N_m =  survivors, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, N_f =  survivors, additive_variance = additive_variance) 
  list(new_state = new_seedbank, survivors = survivors)
}

