#functions for the dynamic programing of the non-spatial model of herb-resistance 
## A crucial part of the dynamic program is the ability to collapse and build normal distributions to and from a given state
## we also want to be confident that a given distirbution is normal
library(lhs)
library(gam)
library(gbm)
library(class)
library(colorspace)
# find mean and sd from a given distirbution
get_mean_sd <- function(dist, eval_points, dg){
  total_sum = sum(dist) * dg
  approx_mean = sum((dist / total_sum) * eval_points * dg)
  approx_sd = sqrt(sum(((eval_points - approx_mean)^2) * (dist / total_sum) * dg))
  return(list(approx_mean = approx_mean, approx_sd = approx_sd, total_pop = total_sum))
}
#build state from a distirbution mean and sd for 2 level seed bank
build_state_2level <- function(dist_top, dist_bottom, eval_points, dg){
  summary_top = get_mean_sd(dist = dist_top, eval_points = eval_points, dg = dg)
  summary_bottom = get_mean_sd(dist = dist_bottom, eval_points = eval_points, dg = dg)
  state = list(mean_top = summary_top$approx_mean, sd_top = summary_top$approx_sd, pop_top = summary_top$total_pop, 
    mean_bottom = summary_bottom$approx_mean, sd_bottom = summary_bottom$approx_sd, pop_bottom = summary_bottom$total_pop)
  return(state)
}
#build a population distirbution from state for 2 level seed bank
build_pops_2level <- function(state, eval_points){
  top = state$pop_top * dnorm(eval_points, mean = state$mean_top, sd = state$sd_top)
  bottom = state$pop_bottom * dnorm(eval_points, mean = state$mean_bottom, sd = state$sd_bottom)
  return(list(top = top, bottom = bottom))
}
#build state from a distirbution mean and sd for 1 level seed bank
build_state_1level <- function(dist_top, eval_points, dg){
  summary_top = get_mean_sd(dist = dist_top, eval_points = eval_points, dg = dg)
  state = list(mean_top = summary_top$approx_mean, sd_top = summary_top$approx_sd, pop_top = summary_top$total_pop)
  return(state)
}
#build a population distirbution from state for 1 level seed bank
build_pops_1level <- function(state, eval_points){
  top = state$pop_top * dnorm(eval_points, mean = state$mean_top, sd = state$sd_top)
  return(top)
}
#test if a given distirbution is normal 
normality_test <- function(dist, eval_points, dg, ...){
  dist_summary = get_mean_sd(dist = dist, eval_points = eval_points, dg = dg)
  normal_dist_approx = dist_summary$total_pop * dnorm(eval_points, mean = dist_summary$approx_mean, sd = dist_summary$approx_sd)
  qqplot(normal_dist_approx, dist, ...)
  abline(0, 1)
}
#get value of action/state pair with resistance 
get_reward_resist <- function(mean_resistance, time_step, time_horizon, budget, num_surviors_pre_mech, income, yeild_loss, actions_nomech, cost_space_nomech, mech, mech_cost0, 
  mech_cost){
  cost = sum(cost_space_nomech[actions_nomech, ]) + mech * (mech_cost0 + mech_cost * num_surviors_pre_mech)
  if(cost > budget){
    return(-Inf)
  }else{
    if(time_step == time_horizon){
      return(-mean_resistance)
    }else{
      return(0)
    }
  }
}
#update the state from a current state to the next one using the proccess model, get next state and reward for current state from the calculation 
#requires a call to fecundity = FECUNDITY_CLOSURE() to create the fecundity() function
state_reward_next_econ <- function(current_state, eval_object, action, sub_action, seed_survival, germination, pro_exposed, max_sur, sur0, sur_cost_resist, 
  herb_effect, survive_resist, fec_max, dense_depend_fec, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, fec_function, dg){
  
  seedbank_current = build_pops_1level(current_state, eval_object$seed) #construct a population from the state
  #calculate the population progression 1 year given the current population and the actions taken.  
  new_seedbank = (seedbank_current - seedbank_current *  sub_action$plow[action['plow']]) * seed_survival 
  new_plants = new_seedbank * germination * sub_action$crop_sur[action['crop']]   
  seedbank_post_germ = new_seedbank * (1 - germination)
  survivors_herb = pro_exposed * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground - 
    sub_action$herb[action['herb']] * (herb_effect - pmin(herb_effect, survive_resist * eval_object$above_ground))))))
  survivors_noherb = (1 - pro_exposed) * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground))))
  survivors_joint = survivors_herb + survivors_noherb #pre-mech population used to calculate cost of mech control
  survivors_joint_mech = survivors_joint * sub_action$mech_sur[action['mech']] #post-mech pop used to caclulate the income from being in that state
  seedbank_next = seedbank_post_germ + fec_function(N_m = survivors_joint_mech, fec_max = fec_max, dense_depend_fec = dense_depend_fec, crop_effect_fec = sub_action$crop_fec[action['crop']], 
    density_effect_fec = sub_action$dens_fec[action['dens']], dg = dg) 
  #take these populations and new seedbank and turn them into rewards and states for the value function calculations 
  income = max(0, income0[action['crop']] - yeild_loss[action['crop']] * sum(survivors_joint_mech) * dg) #money made
  cost_nomech = cost_space_nomech$herb[action['herb']] + cost_space_nomech$crop[action['crop']] + cost_space_nomech$plow[action['plow']] + cost_space_nomech$dens[action['dens']] #non mech cost
  #mech cost
  if(action['mech'] == 1){
    cost_mech = 0
  }else{
    cost_mech = mech_cost0 + mech_cost * sum(survivors_joint) * dg
  }
  
  return(list(next_state = build_state_1level(dist_top = seedbank_next, eval_points = eval_object$seed, dg = dg), reward = income - cost_nomech - cost_mech))
}

#update the state from a current state to the next one using the proccess model, getting the reward along the way, for final time step, T, so only interested in the reward, 
#no need to return future state or run the fecundity function
#requires a call to fecundity = FECUNDITY_CLOSURE() to create the fecundity() function
state_reward_next_econ_T <- function(current_state, eval_object, action, sub_action, seed_survival, germination, pro_exposed, max_sur, sur0, sur_cost_resist, 
  herb_effect, survive_resist, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, dg){
  
  seedbank_current = build_pops_1level(current_state, eval_object$seed) #construct a population from the state
  #calculate the population progression 1 year given the current population and the actions taken.  
  new_seedbank = (seedbank_current - seedbank_current *  sub_action$plow[action['plow']]) * seed_survival 
  new_plants = new_seedbank * germination * sub_action$crop_sur[action['crop']]   
  seedbank_post_germ = new_seedbank * (1 - germination)
  survivors_herb = pro_exposed * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground - 
    sub_action$herb[action['herb']] * (herb_effect - pmin(herb_effect, survive_resist * eval_object$above_ground))))))
  survivors_noherb = (1 - pro_exposed) * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground))))
  survivors_joint = survivors_herb + survivors_noherb #pre-mech population used to calculate cost of mech control
  survivors_joint_mech = survivors_joint * sub_action$mech_sur[action['mech']] #post-mech pop used to caclulate the income from being in that state
  #take these populations and calculate the reward 
  income = max(0, income0[action['crop']] - yeild_loss[action['crop']] * sum(survivors_joint_mech) * dg) #money made
  cost_nomech = cost_space_nomech$herb[action['herb']] + cost_space_nomech$crop[action['crop']] + cost_space_nomech$plow[action['plow']] + cost_space_nomech$dens[action['dens']] #non mech cost
  #mech cost
  if(action['mech'] == 1){
    cost_mech = 0
  }else{
    cost_mech = mech_cost0 + mech_cost * sum(survivors_joint) * dg
  }
  
  return(income - cost_nomech - cost_mech)
}
#For every smapled state I need to find the value of best action return that value, for final time step T
sampled_state_2_value_T <- function(current_state, eval_object, action_space, sub_action, seed_survival, germination, pro_exposed, max_sur, sur0, sur_cost_resist, 
  herb_effect, survive_resist, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, dg){
  
  best_action  = 0
  current_max_value = -Inf
  #get the value for each action a, given the current state
  for(a in 1:dim(action_space)[1]){
    reward_a = state_reward_next_econ_T(current_state = current_state, eval_object = eval_object, action = action_space[a, ], sub_action = sub_action, 
      seed_survival = seed_survival, germination = germination, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, 
      sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, 
      mech_cost0 = mech_cost0, mech_cost = mech_cost, dg = dg)
    
    if(reward_a > current_max_value){
      current_max_value = reward_a
      best_action = a
    }
  }
   
  return(list(value = current_max_value, best_action = best_action))
}
  
#For every smapled state I need to find the value of best action return that value, for a non-final timestep
sampled_state_2_value_gam <- function(current_state, value_surface_t1, discount_factor, eval_object, action_space, sub_action, seed_survival, germination, pro_exposed, max_sur, 
  sur0, sur_cost_resist, herb_effect, survive_resist, fec_max, dense_depend_fec, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, fec_function, dg){
  
  best_action = 0
  current_max_value = -Inf
  #get the value for each action a, given the current state
  for(a in 1:dim(action_space)[1]){
    state_reward = state_reward_next_econ(current_state = current_state, eval_object = eval_object, action = action_space[a, ], sub_action = sub_action, 
      seed_survival = seed_survival, germination = germination, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, 
      herb_effect = herb_effect, survive_resist = survive_resist, fec_max = fec_max, dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yeild_loss, 
      cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, fec_function = fec_function, dg = dg)
   
    value_a = state_reward$reward + discount_factor * predict(value_surface_t1, data.frame(mean_g = state_reward$next_state$mean_top, pop = state_reward$next_state$pop_top)) 
    
    if(value_a > current_max_value){
      current_max_value = value_a
      best_action = a
    }
  }
   
  return(list(value = current_max_value, best_action = best_action))
}

sampled_state_2_value_BRT <- function(current_state, value_surface_t1, discount_factor, eval_object, action_space, sub_action, seed_survival, germination, pro_exposed, max_sur, 
  sur0, sur_cost_resist, herb_effect, survive_resist, fec_max, dense_depend_fec, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, fec_function, dg){
  
  best_action = 0
  current_max_value = -Inf
  #get the value for each action a, given the current state
  for(a in 1:dim(action_space)[1]){
    state_reward = state_reward_next_econ(current_state = current_state, eval_object = eval_object, action = action_space[a, ], sub_action = sub_action, 
      seed_survival = seed_survival, germination = germination, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, 
      herb_effect = herb_effect, survive_resist = survive_resist, fec_max = fec_max, dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yeild_loss, 
      cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, fec_function = fec_function, dg = dg)
   
    value_a = state_reward$reward + discount_factor * predict(value_surface_t1, newdata = data.frame(mean_g = state_reward$next_state$mean_top, pop = state_reward$next_state$pop_top), 
      n.trees = value_surface_t1$n.trees) 
    
    if(value_a > current_max_value){
      current_max_value = value_a
      best_action = a
    }
  }
   
  return(list(value = current_max_value, best_action = best_action))
} 

#string all these functions together to solve the intergrated weed managment problem 
IMW_dynamic_program <- function(inital_state, germination, seed_survival, eval_object_int, pro_exposed, sur0, sur_cost_resist, effect_herb, survive_resist, max_sur, fec_max, fec0, 
  fec_cost, offspring_sd, dense_depend_fec, dg, burnin, sub_action, action_space, time_horizon, discount_factor, income0, yeild_loss, cost_space_nomech, mech_cost0, 
  mech_cost, num_samples, burnin_test_out_loc, burnin_test_out_name = 'burnin_test_output.pdf', BRT_para_list = NA){
    #1. run the model for a burn in period to get the population at some kind of equlibrium under no herbicide, also use this as estimate of lower g that can be expected, 
    noherb_action_seq = cbind(herb = rep(1, burnin), crop = rep(1, burnin), mech = rep(1, burnin), dens = rep(1, burnin), plow = rep(1, burnin))	  
    inital_pop_noherb = multi_iteration(seedbank_initial = inital_state, germination = germination, seed_survival = seed_survival, eval_object = eval_object_int, pro_exposed = pro_exposed,
      sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, 
      offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, dg = dg, num_iter = burnin, sub_action = sub_action, action_seq = noherb_action_seq)
    
    #test the sd stabalises and get the number of individuals over time to estimate the upper population limit.  
    noherb_burnin_list = apply(inital_pop_noherb, MARGIN = 1, FUN = function(x) get_mean_sd(dist = x, eval_points = eval_object$seed, dg = dg))
    #get values over time to test the means, SDs and populations have stablised
    noherb_mean = sapply(noherb_burnin_list, FUN = function(x) x$approx_mean)
    noherb_sd = sapply(noherb_burnin_list, FUN = function(x) x$approx_sd)
    noherb_pop = sapply(noherb_burnin_list, FUN = function(x) x$total_pop)
    #3. run the model for a second burn in period with herbicide applied every year, use this to get upper limit of g for state cost_space_nomech
    herb_action_seq = cbind(herb = rep(2, burnin), crop = rep(1, burnin), mech = rep(1, burnin), dens = rep(1, burnin), plow = rep(1, burnin))	  
    #get values over time to test mean, SD and populations have stabilised 
    inital_pop_herb = multi_iteration(seedbank_initial = inital_state, germination = germination, seed_survival = seed_survival, eval_object = eval_object_int, pro_exposed = pro_exposed,
      sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, 
      offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, dg = dg, num_iter = burnin, sub_action = sub_action, action_seq = herb_action_seq)
    
    herb_burnin_list = apply(inital_pop_herb, MARGIN = 1, FUN = function(x) get_mean_sd(dist = x, eval_points = eval_object$seed, dg = dg))
    herb_mean = sapply(herb_burnin_list, FUN = function(x) x$approx_mean)
    herb_sd = sapply(herb_burnin_list, FUN = function(x) x$approx_sd)
    herb_pop = sapply(herb_burnin_list, FUN = function(x) x$total_pop)
    #make a dataframe fro printing
    print('### values from burn in period under herbicid and no herbicide ###')
    print(cbind(mean_herb = herb_mean, mean_noherb = noherb_mean, sd_herb = herb_sd, sd_noherb = noherb_sd, pop_herb = herb_pop, pop_noherb = noherb_pop))
    
    intial_loc = getwd()
    setwd(burnin_test_out_loc)
    pdf(file = burnin_test_out_name, width = 5, height = 15)#open up plotting window to save the diagnostic plots
      max_mean = max(c(herb_mean, noherb_mean))
      min_mean = min(c(herb_mean, noherb_mean))
      max_sd = max(c(herb_sd, noherb_sd))
      min_sd = min(c(herb_sd, noherb_sd))
      max_pop = max(c(herb_pop, noherb_pop))
      
      par(mfcol = c(3, 1))
      #plot mea over time
      plot(x = seq_along(herb_mean), y = herb_mean, bty = 'n', type = 'l', xlab = 'time', ylab = 'mean g', ylim = c(min_mean, max_mean), col = 'red')
      points(x = seq_along(herb_mean), y = herb_mean, col = 'red', pch = 19)
      points(x = seq_along(noherb_mean), y = noherb_mean, pch = 19, col = 'blue')
      lines(x = seq_along(noherb_mean), y = noherb_mean, col = 'blue')
      #plot sd over time
      plot(x = seq_along(herb_sd), y = herb_sd, bty = 'n', type = 'l', xlab = 'time', ylab = 'sd g', ylim = c(min_sd, max_sd), col = 'red')
      points(x = seq_along(herb_sd), y = herb_sd, col = 'red', pch = 19)
      points(x = seq_along(noherb_sd), y = noherb_sd, pch = 19, col = 'blue')
      lines(x = seq_along(noherb_sd), y = noherb_sd, col = 'blue')
      #plot population over time
      plot(x = seq_along(herb_pop), y = herb_pop, bty = 'n', type = 'l', xlab = 'time', ylab = 'sd g', ylim = c(0, max_pop), col = 'red')
      points(x = seq_along(herb_pop), y = herb_pop, col = 'red', pch = 19)
      points(x = seq_along(noherb_pop), y = noherb_pop, pch = 19, col = 'blue')
      lines(x = seq_along(noherb_pop), y = noherb_pop, col = 'blue')
      legend(x = 'bottomright', legend = c('herb', 'noherb'), col = c('red', 'blue'), lwd = 1, bty = 'n')
      #do the qqplots to test for normality
      for(i in 1:dim(inital_pop_herb)[1]){
	par(mfcol = c(2, 1))
	normality_test(dist = inital_pop_herb[i, ], eval_points = eval_object$seed, dg = dg, main = paste0('herb | time: ', i))
	normality_test(dist = inital_pop_noherb[i, ], eval_points = eval_object$seed, dg = dg, main = paste0('no_herb | time: ', i))
      }
    dev.off()
    setwd(intial_loc)
    #use these end points to set up the state space and re-define the evaluation points to make the domian fit the limits of the state space more closley 
    sd_est = mean(c(herb_sd[burnin], noherb_sd[burnin]))
    mean_eval_lower = floor(noherb_mean[burnin] - 9 * sd_est)
    mean_eval_upper = ceiling(herb_mean[burnin] + 9 * sd_est)
    
    mean_state_lower = floor(noherb_mean[burnin] - 3 * sd_est)
    mean_state_upper = ceiling(herb_mean[burnin] + 3 * sd_est)
    pop_state_upper = max(c(herb_pop, noherb_pop))
 
    #rebuild the evaluation points object to be smaller (possibly), set it as mean value +- 3*sd 
    eval_obj_mod = eval_points_builder(lower_eval_point = mean_eval_lower, upper_eval_point = mean_eval_upper, resolution = dg, seed_expantion = ceiling(4 * offspring_sd)) 
    sample_points = improvedLHS(n = num_samples, k = 2, dup = 3)
    sample_points[, 1] = sample_points[, 1] * (mean_state_upper - mean_state_lower) + mean_state_lower
    sample_points[, 2] = sample_points[, 2] * pop_state_upper
    sample_points = rbind(sample_points, c(mean_state_lower, 0.5), c(mean_state_lower, pop_state_upper), c(mean_state_upper, 0.5), c(mean_state_upper, pop_state_upper)) #adds the 4 extream corners of the state space to make sure the gam does not end up extrapolating outside the range of the data
    sample_points_df = data.frame(mean_g = sample_points[, 1], pop = sample_points[, 2])#turn the sample points to a data frame since gam() likes dataframes the most
    #create the Q object to hold the action state values, along with value surface and best actions at sample points to help approximate the policy 
    Q = list(parameters = list(params = c(germination = germination, seed_survival = seed_survival, pro_exposed = pro_exposed, sur0 = sur0, sur_cost_resist = sur_cost_resist, 
      effect_herb = effect_herb, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, pop_sd = sd_est, 
      dense_depend_fec = dense_depend_fec, time_horizon = time_horizon, discount_factor = discount_factor, mech_cost0 = mech_cost0, mech_cost = mech_cost, dg = dg), 
      eval_obj = eval_obj_mod, sub_action = sub_action, action_space = action_space, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech), 
      DP = list()) 
    #find the value of each sampled point in the sate space, given best action taken, also what the best action was
    value_action_sample_points = apply(sample_points, MARGIN = 1, FUN = function(x){
      sampled_state_2_value_T(current_state = list(mean_top = x[1], sd_top = sd_est, pop_top = x[2]), eval_object = eval_obj_mod, action_space = action_space, sub_action = sub_action,
	seed_survival = seed_survival, germination = germination, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, 
	herb_effect = effect_herb, survive_resist = survive_resist, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, 
	mech_cost = mech_cost, dg = dg)
    })
    
    #add the values to the sampled points in state space
    sample_points_df$values = sapply(value_action_sample_points, FUN = function(x) x$value)
    sample_points_df$policy = sapply(value_action_sample_points, FUN = function(x) x$best_action)
    
    #check if BRT parameters set, or if blank, in which case fit gam 
    if(length(BRT_para_list) == 1){
      cat('calculating time step: ', 1, '\n')
      value_surface_t = suppressWarnings(gam(values ~ s(mean_g, spar = 0.01)  + s(pop, spar = 0.01), data = sample_points_df))
      Q$DP[[1]] = list(sample_points_df = sample_points_df, value_surface = value_surface_t)
      #make the fecundity function on the new evaluation points 
      fecundity = fecundity_closure(eval_points_object = eval_obj_mod, offspring_sd = offspring_sd, fec0 = fec0, fec_cost = fec_cost)
    
      #start stepping back in time 
      for(i in 2:time_horizon){
	cat('calculating time step: ', i, '\n')
	value_action_sample_points = apply(sample_points, MARGIN = 1, FUN = function(x){
	  sampled_state_2_value_gam(current_state = list(mean_top = x[1], sd_top = sd_est, pop_top = x[2]), value_surface_t1 = Q$DP[[i - 1]]$value_surface, discount_factor = discount_factor, 
	    eval_object = eval_obj_mod, action_space = action_space, sub_action = sub_action, seed_survival = seed_survival, germination = germination, pro_exposed = pro_exposed, 
	    max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, fec_max = fec_max, 
	    dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, 
	    fec_function = fecundity, dg = dg)
	})
	sample_points_df$values = sapply(value_action_sample_points, FUN = function(x) x$value)
	sample_points_df$policy = sapply(value_action_sample_points, FUN = function(x) x$best_action)
	value_surface_t = suppressWarnings(gam(values ~ s(mean_g, spar = 0.01) + s(pop, spar = 0.01), data = sample_points_df))
	Q$DP[[i]] = list(sample_points_df = sample_points_df, value_surface = value_surface_t)
      }
      return(Q)
    }else{
      cat('calculating time step: ', 1, '\n')
      value_surface_t = gbm.fit(x = cbind(sample_points_df$mean_g, sample_points_df$pop), y = sample_points_df$values, distribution = 'gaussian', n.trees = BRT_para_list$num_trees, interaction.depth = BRT_para_list$int_depth, 
	shrinkage = BRT_para_list$learn_rate, bag.fraction = 1, keep.data = FALSE, verbose = FALSE)
      Q$DP[[1]] = list(sample_points_df = sample_points_df, value_surface = value_surface_t)
      #make the fecundity function on the new evaluation points 
      fecundity = fecundity_closure(eval_points_object = eval_obj_mod, offspring_sd = offspring_sd, fec0 = fec0, fec_cost = fec_cost)
      for(i in 2:time_horizon){
	cat('calculating time step: ', i, '\n')
	value_action_sample_points = apply(sample_points, MARGIN = 1, FUN = function(x){
	  sampled_state_2_value_BRT(current_state = list(mean_top = x[1], sd_top = sd_est, pop_top = x[2]), value_surface_t1 = Q$DP[[i - 1]]$value_surface, discount_factor = discount_factor, 
	    eval_object = eval_obj_mod, action_space = action_space, sub_action = sub_action, seed_survival = seed_survival, germination = germination, pro_exposed = pro_exposed, 
	    max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, fec_max = fec_max, 
	    dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, 
	    fec_function = fecundity, dg = dg)
	})
	sample_points_df$values = sapply(value_action_sample_points, FUN = function(x) x$value)
	sample_points_df$policy = sapply(value_action_sample_points, FUN = function(x) x$best_action)
	value_surface_t = gbm.fit(x = cbind(sample_points_df$mean_g, sample_points_df$pop), y = sample_points_df$values, distribution = 'gaussian', n.trees = BRT_para_list$num_trees, interaction.depth = BRT_para_list$int_depth, 
	  shrinkage = BRT_para_list$learn_rate, bag.fraction = 1, keep.data = FALSE, verbose = FALSE)
	Q$DP[[i]] = list(sample_points_df = sample_points_df, value_surface = value_surface_t)
      }
      return(Q)
    }
}
 
#simulate a policy from a set of starting points, record the state every time step 
policy_simulator <- function(policy, model_pars, time_period, start_points){
  mean_g_range = range(policy$sample_points_df$mean_g)
  pop_range = range(policy$sample_points_df$pop)
  #set up the scales training points for knn interpilation of the policy over state space
  mean_g_scaled = (policy$sample_points_df$mean_g - mean_g_range[1]) / (mean_g_range[2] - mean_g_range[1])
  pop_scaled = (policy$sample_points_df$pop - pop_range[1]) / (pop_range[2] - pop_range[1])
  #create fecundity function 
  fecundity = fecundity_closure(eval_points_object = model_pars$eval_obj, offspring_sd = model_pars$params['offspring_sd'], fec0 = model_pars$params['fec0'], 
    fec_cost = model_pars$params['fec_cost'])
  #for every pointed sampled in the state space find the fate of the population at the time horizon for following that 
  
  results = apply(start_points, MARGIN = 1, FUN = function(s){
  
    #objects to fill and return
    action_seq = matrix(NA, nrow = time_period, ncol = dim(model_pars$action_space)[2] + 1)
    colnames(action_seq) <- c('act_num', colnames(model_pars$action_space))
    state_over_time = list()
   
    #build the action vector
    new_point = data.frame(mean_g = (s['mean_g'] - mean_g_range[1]) / (mean_g_range[2] - mean_g_range[1]),
      pop = (s['pop'] - pop_range[1]) / (pop_range[2] - pop_range[1]))
    action_seq[1, 'act_num'] = as.numeric(as.character(knn1(train = data.frame(mean_g = mean_g_scaled, pop = pop_scaled), test = new_point, cl = policy$sample_points_df$policy)))
    action_seq[1, 2:(dim(model_pars$action_space)[2] + 1)] = model_pars$action_space[action_seq[1, 'act_num'], ]
    #set the intial population
    state_over_time[[1]] = state_reward_next_econ(current_state = list(mean_top = s['mean_g'], sd_top = model_pars$params['pop_sd'], pop_top = s['pop']), 
      eval_object = model_pars$eval_obj, action = model_pars$action_space[action_seq[1, 'act_num'], ], sub_action = model_pars$sub_action, seed_survival = model_pars$params['seed_survival'], 
      germination = model_pars$params['germination'], pro_exposed = model_pars$params['pro_exposed'], max_sur = model_pars$params['max_sur'], sur0 = model_pars$params['sur0'], 
      sur_cost_resist = model_pars$params['sur_cost_resist'], herb_effect = model_pars$params['effect_herb'], survive_resist = model_pars$params['survive_resist'], 
      fec_max = model_pars$params['fec_max'], dense_depend_fec = model_pars$params['dense_depend_fec'], income0 = model_pars$income0, yeild_loss = model_pars$yeild_loss, 
      cost_space_nomech = model_pars$cost_space_nomech, mech_cost0 = model_pars$params['mech_cost0'], mech_cost = model_pars$params['mech_cost'], fec_function = fecundity, 
      dg = model_pars$params['dg'])
    #simulate over timesteps
    for(ts in 2:time_period){
      #test to make sure the new population stays in the state space, break the loop
      if(state_over_time[[ts - 1]]$next_state$mean_top > mean_g_range[2] | state_over_time[[ts - 1]]$next_state$mean_top < mean_g_range[1] | 
	state_over_time[[ts - 1]]$next_state$pop_top > pop_range[2] | state_over_time[[ts - 1]]$next_state$pop_top < pop_range[1]) break
      #scale the new state so that knn1 works
      new_point = data.frame(mean_g = (state_over_time[[ts - 1]]$next_state$mean_top - mean_g_range[1]) / (mean_g_range[2] - mean_g_range[1]), 
	pop = (state_over_time[[ts - 1]]$next_state$pop_top - pop_range[1]) / (pop_range[2] - pop_range[1]))
      #find the best action by interpolation of policy using knn1
      action_seq[ts, 'act_num'] = as.numeric(as.character(knn1(train = data.frame(mean_g = mean_g_scaled, pop = pop_scaled), test = new_point, cl = policy$sample_points_df$policy)))
      action_seq[ts, 2:(dim(model_pars$action_space)[2] + 1)] = model_pars$action_space[action_seq[ts, 'act_num'], ]
      
      state_over_time[[ts]] = state_reward_next_econ(current_state = state_over_time[[ts - 1]]$next_state, eval_object = model_pars$eval_obj, 
	action = model_pars$action_space[action_seq[ts, 'act_num'], ], sub_action = model_pars$sub_action, seed_survival = model_pars$params['seed_survival'], 
	germination = model_pars$params['germination'], pro_exposed = model_pars$params['pro_exposed'], max_sur = model_pars$params['max_sur'], sur0 = model_pars$params['sur0'], 
	sur_cost_resist = model_pars$params['sur_cost_resist'], herb_effect = model_pars$params['effect_herb'], survive_resist = model_pars$params['survive_resist'], 
	fec_max = model_pars$params['fec_max'], dense_depend_fec = model_pars$params['dense_depend_fec'], income0 = model_pars$income0, yeild_loss = model_pars$yeild_loss, 
	cost_space_nomech = model_pars$cost_space_nomech, mech_cost0 = model_pars$params['mech_cost0'], mech_cost = model_pars$params['mech_cost'], fec_function = fecundity, 
	dg = model_pars$params['dg'])
    }
  
    return(list(start_point = s, action_seq = action_seq, state_over_time = state_over_time))
  })  
  
  return(results)
}


