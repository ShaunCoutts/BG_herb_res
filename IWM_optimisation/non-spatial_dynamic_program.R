#functions for the dynamic programing of the non-spatial model of herb-resistance 
## A crucial part of the dynamic program is the ability to collapse and build normal distributions to and from a given state
## we also want to be confident that a given distirbution is normal
# find mean and sd from a given distirbution
get_mean_sd <- function(dist, eval_points, dg){
  total_sum = sum(dist * dg)
  approx_mean = sum((dist / total_sum) * eval_points * dg)
  approx_sd = sqrt(sum(((eval_points - approx_mean)^2) * (dist / total_sum) * dg))
  return(list(approx_mean = approx_mean, approx_sd = approx_sd, total_pop = total_sum))
}
#build state from a distirbution mean and sd
build_state <- function(dist_top, dist_bottom, eval_points, dg){
    summary_top = get_mean_sd(dist = dist_top, eval_points = eval_points, dg = dg)
    summary_bottom = get_mean_sd(dist = dist_bottom, eval_points = eval_points, dg = dg)
    state = list(mean_top = summary_top$approx_mean, sd_top = summary_top$approx_sd, pop_top = summary_top$total_pop, 
      mean_bottom = summary_bottom$approx_mean, sd_bottom = summary_bottom$approx_sd, pop_bottom = summary_bottom$total_pop)
    return(state)
}
#build a population distirbution from state
build_pops <- function(state, eval_points){
  top = state$pop_top * dnorm(eval_points, mean = state$mean_top, sd = state$sd_top)
  bottom = state$pop_bottom * dnorm(eval_points, mean = state$mean_bottom, sd = state$sd_bottom)
  return(list(top = top, bottom = bottom))
}
#test if a given distirbution is normal 
normality_test <- function(dist, eval_points, dg){
    dist_summary = get_mean_sd(dist = dist, eval_points = eval_points, dg = dg)
    normal_dist_approx = dist_summary$total_pop * dnorm(eval_points, mean = dist_summary$approx_mean, sd = dist_summary$approx_sd)
    qqplot(normal_dist_approx, dist)
    abline(0, 1)
}

## use of the system model to project from one state to another, also to obtain intial population
#create the intial population from a starting normal distribution in order to get the standard deviation to apply  
create_intial_population <- function(intial_population){

}
















## creat a node with all the information needed.
node <- function(act, parent, alt_profit, seed_survival, base_profit, weed_impact, germination, eval_object, sur0, sur_cost_resist, herb_effect, 
  survive_resist, max_sur, ceiling_pop, fec_max, fec0, fec_cost, additive_variance, density_cutoff, dg, dis_rate){
  if(act == 3){
    new_state = seedbank_update(parent$state, seed_survival)
    value = parent$value + ((base_profit - alt_profit) / ((1 + dis_rate) ^ (parent$time_step + 1)))
  }else{
    N1 = state_transition(seedbank_current = parent$state, germination = germination, eval_object = eval_object, herb_rate = A_herb[act], sur0 = sur0, sur_cost_resist = sur_cost_resist, 
      herb_effect = herb_effect, survive_resist = sur_cost_resist, max_sur = max_sur, ceiling_pop = ceiling_pop, seed_survival = seed_survival, fec_max = fec_max, fec0 = fec0, 
      fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff)
    new_state = N1$new_state
    value = parent$value + ((base_profit - reward(act, N1$survivors, base_profit, weed_impact, dg)) / ((1 + dis_rate) ^ (parent$time_step + 1)))
  }
  
  list(acts = c(act, parent$acts), time_step = parent$time_step + 1, value = value, state = new_state)
}

intial_node <- function(state0){
  list(acts = 0, time_step = 0, value = 0, state = state0)
}

add_node <- function(node_list, act, parent_index, alt_profit, seed_survival, base_profit, weed_impact, germination, eval_object, sur0, sur_cost_resist, herb_effect, 
  survive_resist, max_sur, ceiling_pop, fec_max, fec0, fec_cost, additive_variance, density_cutoff, dg, dis_rate){
  end = length(node_list)
  node_list[[end + 1]] = node(act = act, parent = node_list[[parent_index]], alt_profit = alt_profit, seed_survival = seed_survival, base_profit = base_profit,
    weed_impact = weed_impact, germination = germination, eval_object = eval_object, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist,
    max_sur = max_sur, ceiling_pop = ceiling_pop, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff, 
    dg = dg, dis_rate = dis_rate)
    
    node_list
}

reward <- function(act, survivors, base_profit, weed_impact, dg){
  N = sum(survivors) * dg
  max(0, base_profit - N * weed_impact) - A_cost[act] 
}

seedbank_update <- function(seedbank_current, seed_survival){
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


simulate_solution <- function(solution_list, eval_object, pause = 0.5, ...){
  max_value = max(sapply(solution_list, FUN = function(x) x$state))
  actions = c('No Herbicide', 'Herbicide', 'Alt Crop')
  cols = c('skyblue', 'red', 'green')
  for(i in 1:length(solution_list)){
    title = paste0('Action: ', ifelse(solution_list[[i]]$acts[1] == 0, 'intial step', actions[solution_list[[i]]$acts[1]]), 'Value = ', solution_list[[i]]$value, '\ntime =', solution_list[[i]]$time_step)
    plot(eval_object$seed, eval_object$seed, type = 'n', ylim = c(0, max_value), bty = 'n', xlab = 'resistance score', main = title, ...)
    polygon(x = eval_object$seed, y = solution_list[[i]]$state, col = ifelse(solution_list[[i]]$acts[1] == 0, 'white', cols[solution_list[[i]]$acts[1]]))
    Sys.sleep(pause)
  }
}



## some tests to check things are working 
tests <- function(){
  test_out = character()
  #test the mean and sd recovery
  test_mean = 0.28
  test_sd = 3.6
  dx = 0.5
  N = 168.9
  x = seq(-50, 50, dx)
  y = N * dnorm(x, test_mean, test_sd)
  means_sd = get_mean_sd(dist = y, eval_points = x, dg = dx)
  test_out[1] <- ifelse(means_sd$approx_mean == test_mean & means_sd$approx_sd == test_sd & means_sd$total_pop == N, 'GET_MEAN_SD() worked', 'GET_MEAN_SD() failed')

  #make a plot to see if the normality test produces a plot, first make a  comparasion to a non-normal distribution (log-normal)
  test_dist_non_norm = N * dlnorm(x, meanlog = 3.2, sdlog = 1)
  test_dist_norm = N * dnorm(x, mean = 3.2, sd = 1)
  par(mfrow = c(1, 2))
  normality_test(test_dist_non_norm, x, dx)
  normality_test(test_dist_norm, x, dx)
  
  print(test_out)
}
tests()
