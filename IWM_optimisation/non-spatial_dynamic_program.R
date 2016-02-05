#functions for the dynamic programing of the non-spatial model of herb-resistance 
## A crucial part of the dynamic program is the ability to collapse and build normal distributions to and from a given state
## we also want to be confident that a given distirbution is normal
library(lhs)
library(gam)
library(colorspace)
# find mean and sd from a given distirbution
get_mean_sd <- function(dist, eval_points, dg){
  total_sum = sum(dist * dg)
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
state_reward_next_econ <- function(current_state, eval_object, action, sub_action, seed_survival, germination, density_cutoff, pro_exposed, max_sur, sur0, sur_cost_resist, 
  herb_effect, survive_resist, fec_max, fec0, fec_cost, offspring_sd, dense_depend_fec, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, dg){
  
  seedbank_current = build_pops_1level(current_state, eval_object$seed) #construct a population from the state
  #calculate the population progression 1 year given the current population and the actions taken.  
  new_seedbank = (seedbank_current - seedbank_current *  sub_action$plow[action['plow']]) * seed_survival 
  new_plants = new_seedbank * germination * sub_action$crop_sur[action['crop']]   
  seedbank_post_germ = new_seedbank * (1 - germination)
  eval_object = eval_points_update(eval_points_object = eval_object, new_plants, density_cutoff = density_cutoff) #update evaluation window
  survivors_herb = pro_exposed * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground - 
    sub_action$herb[action['herb']] * (herb_effect - pmin(herb_effect, survive_resist * eval_object$above_ground))))))
  survivors_noherb = (1 - pro_exposed) * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground))))
  survivors_joint = survivors_herb + survivors_noherb #pre-mech population used to calculate cost of mech control
  survivors_joint_mech = survivors_joint * sub_action$mech_sur[action['mech']] #post-mech pop used to caclulate the income from being in that state
  seedbank_next = seedbank_post_germ + fecundity(N_m = survivors_joint_mech, eval_points = eval_object$above_ground, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, 
    N_f = survivors_joint_mech, offspring_sd = offspring_sd, seed_eval_points = eval_object$seed, dense_depend_fec = dense_depend_fec, crop_effect_fec = sub_action$crop_fec[action['crop']], 
    density_effect_fec = sub_action$dens_fec[action['dens']], dg = dg) 
  #take these populations and new seedbank and turn them into rewards and states for the value function calculations 
  income = max(0, income0[action['crop']] - yeild_loss[action['crop']] * sum(survivors_joint_mech) * dg) #money made
  cost_nomech = cost_space_nomech$herb[action['herb']] + cost_space_nomech$crop[action['crop']] + cost_space_nomech$plow[action['plow']] + cost_space_nomech$dens[action['dens']] #non mech cost
  #mech cost
  if(action['mech'] == 0){
    cost_mech = 0
  }else{
    cost_mech = mech_cost0 + mech_cost * sum(survivors_joint) * dg
  }
  
  return(list(next_state = build_state_1level(dist_top = seedbank_next, eval_points = eval_object$seed, dg = dg), reward = income - cost_nomech - cost_mech))
}

#update the state from a current state to the next one using the proccess model, getting the reward along the way, for final time step, T, so only interested in the reward, 
#no need to return future state or run the fecundity function
state_reward_next_econ_T <- function(current_state, eval_object, action, sub_action, seed_survival, germination, density_cutoff, pro_exposed, max_sur, sur0, sur_cost_resist, 
  herb_effect, survive_resist, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, dg){
  
  seedbank_current = build_pops_1level(current_state, eval_object$seed) #construct a population from the state
  #calculate the population progression 1 year given the current population and the actions taken.  
  new_seedbank = (seedbank_current - seedbank_current *  sub_action$plow[action['plow']]) * seed_survival 
  new_plants = new_seedbank * germination * sub_action$crop_sur[action['crop']]   
  seedbank_post_germ = new_seedbank * (1 - germination)
  eval_object = eval_points_update(eval_points_object = eval_object, new_plants, density_cutoff = density_cutoff) #update evaluation window
  survivors_herb = pro_exposed * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground - 
    sub_action$herb[action['herb']] * (herb_effect - pmin(herb_effect, survive_resist * eval_object$above_ground))))))
  survivors_noherb = (1 - pro_exposed) * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground))))
  survivors_joint = survivors_herb + survivors_noherb #pre-mech population used to calculate cost of mech control
  survivors_joint_mech = survivors_joint * sub_action$mech_sur[action['mech']] #post-mech pop used to caclulate the income from being in that state
  #take these populations and calculate the reward 
  income = max(0, income0[action['crop']] - yeild_loss[action['crop']] * sum(survivors_joint_mech) * dg) #money made
  cost_nomech = cost_space_nomech$herb[action['herb']] + cost_space_nomech$crop[action['crop']] + cost_space_nomech$plow[action['plow']] + cost_space_nomech$dens[action['dens']] #non mech cost
  #mech cost
  if(action['mech'] == 0){
    cost_mech = 0
  }else{
    cost_mech = mech_cost0 + mech_cost * sum(survivors_joint) * dg
  }
  
  return(income - cost_nomech - cost_mech)
}
#For every smapled state I need to find the value of best action return that value, for final time step T
sampled_state_2_value_T <- function(current_state, eval_object, action_space, sub_action, seed_survival, germination, density_cutoff, pro_exposed, max_sur, sur0, sur_cost_resist, 
  herb_effect, survive_resist, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, dg){
  
  best_action  = 0
  current_max_value = -Inf
  #get the value for each action a, given the current state
  for(a in 1:dim(action_space)[1]){
    reward_a = state_reward_next_econ_T(current_state = current_state, eval_object = eval_object, action = action_space[a, ], sub_action = sub_action, 
      seed_survival = seed_survival, germination = germination, density_cutoff = density_cutoff, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, 
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
sampled_state_2_value <- function(current_state, value_surface_t1, discount_factor, eval_object, action_space, sub_action, seed_survival, germination, density_cutoff, pro_exposed, max_sur, sur0, sur_cost_resist, 
  herb_effect, survive_resist, fec_max, fec0, fec_cost, offspring_sd, dense_depend_fec, income0, yeild_loss, cost_space_nomech, mech_cost0, mech_cost, dg){
  
  best_action = 0
  current_max_value = -Inf
  #get the value for each action a, given the current state
  for(a in 1:dim(action_space)[1]){
    state_reward = state_reward_next_econ(current_state = current_state, eval_object = eval_object, action = action_space[a, ], sub_action = sub_action, 
      seed_survival = seed_survival, germination = germination, density_cutoff = density_cutoff, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, 
      sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, 
      offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, 
      mech_cost0 = mech_cost0, mech_cost = mech_cost, dg = dg)
   
    value_a = state_reward$reward + discount_factor * predict(value_surface_t1, data.frame(mean_g = state_reward$next_state$mean_top, pop = state_reward$next_state$pop_top)) 
    
    if(value_a > current_max_value){
      current_max_value = value_a
      best_action = a
    }
  }
   
  return(list(value = current_max_value, best_action = best_action))
}
  

#string all these functions together to solve the intergrated weed managment problem 
IMW_dynamic_program <- function(inital_state, germination, seed_survival, eval_object_int, pro_exposed, sur0, sur_cost_resist, effect_herb, survive_resist, max_sur, fec_max, fec0, 
  fec_cost, offspring_sd, dense_depend_fec, density_cutoff, dg, burnin, sub_action, action_space, time_horizon, discount_factor, income0, yeild_loss, cost_space_nomech, mech_cost0, 
  mech_cost, burnin_test_out_loc, burnin_test_out_name = 'burnin_test_output.pdf'){
    #1. run the model for a burn in period to get the population at some kind of equlibrium under no herbicide, also use this as estimate of lower g that can be expected, 
    noherb_action_seq = cbind(herb = rep(1, burnin), crop = rep(1, burnin), mech = rep(1, burnin), dens = rep(1, burnin), plow = rep(1, burnin))	  
 
    inital_pop_noherb = multi_iteration(seedbank_initial = inital_state, germination = germination, seed_survival = seed_survival, eval_object = eval_object_int, pro_exposed = pro_exposed,
      sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, 
      offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, density_cutoff = density_cutoff, dg = dg, num_iter = burnin, sub_action = sub_action, action_seq = noherb_action_seq)
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
      offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, density_cutoff = density_cutoff, dg = dg, num_iter = burnin, sub_action = sub_action, action_seq = herb_action_seq)
    
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
    mean_lower = floor(noherb_mean[burnin] - 3 * sd_est)
    mean_upper = ceiling(herb_mean[burnin] + 3 * sd_est)
    pop_upper = max(c(herb_pop, noherb_pop))
    
    #rebuild the evaluation points object to be smaller (possibly), set it as mean value +- 3*sd 
    eval_obj_mod = eval_points_builder(lower_eval_point = mean_lower, upper_eval_point = mean_upper, resolution = dg, seed_expantion = ceiling(3 * offspring_sd)) 
    sample_points = improvedLHS(n = 100, k = 2, dup = 3)
    sample_points[, 1] = sample_points[, 1] * (mean_upper - mean_lower) + mean_lower
    sample_points[, 2] = sample_points[, 2] * pop_upper
    sample_points = rbind(sample_points, c(mean_lower, 0.5), c(mean_lower, pop_upper), c(mean_upper, 0.5), c(mean_upper, pop_upper)) #adds the 4 extream corners of the state space to makesure the gamdoes not end up extrapolating outside the range of the data
    sample_points_df = data.frame(mean_g = sample_points[, 1], pop = sample_points[, 2])#turn the sample points to a data frame since gam() likes dataframes the most
    
    Q = list() #create the Q object to hold the action state values, along with value surface and best actions at sample points to help approximate the policy 
    #find the value of each sampled point in the sate space, given best action taken, also what the best action was
    value_action_sample_points = apply(sample_points, MARGIN = 1, FUN = function(x){
      sampled_state_2_value_T(current_state = list(mean_top = x[1], sd_top = sd_est, pop_top = x[2]), eval_object = eval_obj_mod, action_space = action_space, sub_action = sub_action,
	seed_survival = seed_survival, germination = germination, density_cutoff = density_cutoff, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, 
	herb_effect = effect_herb, survive_resist = survive_resist, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, 
	mech_cost = mech_cost, dg = dg)
    })
    #add the values to the sampled points in state space
    sample_points_df$values = sapply(value_action_sample_points, FUN = function(x) x$value)
    sample_points_df$policy = sapply(value_action_sample_points, FUN = function(x) x$best_action)
    value_surface_t = suppressWarnings(gam(values ~ s(mean_g, spar = 0.01)  + s(pop, spar = 0.01), data = sample_points_df))
    Q[[1]] = list(sample_points_df = sample_points_df, value_surface = value_surface_t)
    #start stepping back in time 
    for(i in 2:time_horizon){
      value_action_sample_points = apply(sample_points, MARGIN = 1, FUN = function(x){
	sampled_state_2_value(current_state = list(mean_top = x[1], sd_top = sd_est, pop_top = x[2]), value_surface_t1 = Q[[i - 1]]$value_surface, discount_factor = discount_factor, 
	  eval_object = eval_obj_mod, action_space = action_space, sub_action = sub_action, seed_survival = seed_survival, germination = germination, density_cutoff = density_cutoff, 
	  pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = effect_herb, survive_resist = survive_resist, fec_max = fec_max, 
	  fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yeild_loss, 
	  cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, dg = dg)
      })
      sample_points_df$values = sapply(value_action_sample_points, FUN = function(x) x$value)
      sample_points_df$policy = sapply(value_action_sample_points, FUN = function(x) x$best_action)
      value_surface_t = suppressWarnings(gam(values ~ s(mean_g, spar = 0.01) + s(pop, spar = 0.01), data = sample_points_df))
      Q[[i]] = list(sample_points_df = sample_points_df, value_surface = value_surface_t)
    }
    return(Q)
 }
 
#simple plot of best action in action space on the sampled state space
plot_policy <- function(Q, output_loc, output_name = 'policy_output.pdf'){
  inital_loc = getwd()
  setwd(output_loc)
  col_palett = rainbow_hcl(n = 25, start = 0, end = 200) 
  #set up an evaluation grid to plot the value surface from the gam 
  upper_mean = max(Q[[1]]$sample_points_df$mean_g)
  lower_mean = min(Q[[1]]$sample_points_df$mean_g)
  upper_pop = max(Q[[1]]$sample_points_df$pop)
  lower_pop = min(Q[[1]]$sample_points_df$pop)
  eval_res = 100
  value_grid = data.frame(mean_g = rep.int(seq(lower_mean, upper_mean, length.out = eval_res), eval_res), pop = rep(seq(lower_pop, upper_pop, length.out = eval_res), each = eval_res))#make every combination of evaluation points
  pdf(file = output_name, width = 15, height = 5)
    for(i in seq_along(Q)){
      par(mfrow = c(1, 3))
      plot(x = Q[[i]]$sample_points_df$mean_g, y = Q[[i]]$sample_points_df$pop, type = 'n', xlab = 'mean_g', ylab = 'population', main = paste0('Policy t = T - ', (i - 1)))
      text(x = Q[[i]]$sample_points_df$mean_g , y = Q[[i]]$sample_points_df$pop , labels = paste0('a', Q[[i]]$sample_points_df$policy), 
	col = col_palett[Q[[i]]$sample_points_df$policy])
      #plot the value surface
      image(matrix(predict(Q[[i]]$value_surface, value_grid), ncol = eval_res, byrow = FALSE), xaxt = 'n', xlab = 'mean_g', yaxt = 'n', ylab = 'population')
      axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(value_grid$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1))
      axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(value_grid$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 1))
      #plot fitted VS predicted for spline
      plot(Q[[i]]$sample_points_df$values, predict(Q[[i]]$value_surface, Q[[i]]$sample_points_df), xlab = 'observed', ylab = 'predicted') 
      abline(0, 1)
    }
  dev.off()
  setwd(inital_loc)
}
#run the functions
DP_policy = IMW_dynamic_program(inital_state = inital_state, germination = germination, seed_survival = seed_survival, eval_object_int = eval_object_int, pro_exposed = pro_exposed,
  sur0 = sur0, sur_cost_resist = sur_cost_resist, effect_herb = effect_herb, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost,
  offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, density_cutoff = density_cutoff, dg = dg, burnin = 50, sub_action = sub_action, action_space = action_space, 
  time_horizon = 10, discount_factor = 0.96, income0 = income0, yeild_loss = yeild_loss, cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, 
  burnin_test_out_loc = output_loc)

plot_policy(Q = DP_policy, output_loc = output_loc)

#play around with a few spline fitting tools here to see how they work
# data(mtcars)
# mtcars_norep = mtcars[!duplicated(mtcars$mpg), ]
# mtcars_norep = mtcars_norep[!duplicated(mtcars_norep$qsec), ]
# mtcars_norep = mtcars_norep[!duplicated(mtcars_norep$wt), ]
# #set up some plots to see how well the 
# par(mfrow = c(2,2))
# interp = gam(mpg ~ s(wt, spar = 1) + s(qsec, spar = 1), data = mtcars)
# plot(mtcars$mpg, predict(interp, mtcars), main = 's-spline, spar = 1')
# abline(0, 1)
# 
# interp = gam(mpg ~ s(wt, spar = 0.01) + s(qsec, spar = 0.01), data = mtcars)
# plot(mtcars$mpg, predict(interp, mtcars), main = 's-spline, spar = 0.01')
# abline(0, 1)
# 
# interp = gam(mpg ~ poly(wt, degree = 3, raw = TRUE) + poly(qsec, degree = 3, raw = TRUE), data = mtcars)
# plot(mtcars$mpg, predict(interp, mtcars), main = 'poly-spline, degree = 3')
# abline(0, 1)
# 
# interp = gam(mpg ~ poly(wt, degree = length(mtcars$mpg) - 1, raw = TRUE) + poly(qsec, degree = length(mtcars$mpg) - 1, raw = TRUE), data = mtcars)
# plot(mtcars$mpg, predict(interp, mtcars), main = 's-spline, df = n - 1')
# abline(0, 1)
# 
# 
# #some speed testing to see if the LHS sampling speed changes much between methods
# library(microbenchmark)
# microbenchmark(randomLHS(n = 100, k = 2), improvedLHS(100, k = 2), optimumLHS(n = 100, k = 2), times = 100)
# par(mfrow = c(1, 3))
# plot(randomLHS(n = 100, k = 2))
# plot(improvedLHS(n = 100, k = 2))
# plot(optimumLHS(n = 100, k = 2))
