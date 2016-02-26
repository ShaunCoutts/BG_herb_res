
##Bunch of testing and plotting functions for dynamic program

## SINGLE_ITERATION_PLOT()
##just like SINGLE_INTERATION_1LEVEL(), but plots all the internal distributions for diagnostic purposes
single_iteration_plot <- function(seedbank_current, germination, mech_control, crop_effect_sur, seed_survival, seed_movement, eval_object, pro_exposed, herb_rate, sur0, 
    sur_cost_resist, herb_effect, survive_resist, max_sur, fec_max, dense_depend_fec, crop_effect_fec, density_effect_fec, fec_function, iteration, dg){
  
  new_seedbank = (seedbank_current - seedbank_current *  seed_movement) * seed_survival 
  new_plants = new_seedbank * germination * mech_control * crop_effect_sur   
  seedbank_post_germ = new_seedbank * (1 - germination)
  survivors_herb = pro_exposed * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground - herb_rate * (herb_effect - 
    pmin(herb_effect, survive_resist * eval_object$above_ground))))))
  survivors_noherb = (1 - pro_exposed) * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground))))
  survivors_joint = survivors_herb + survivors_noherb
  new_seeds = fec_function(N_m = survivors_joint, fec_max = fec_max, dense_depend_fec = dense_depend_fec, 
    density_effect_fec = density_effect_fec, crop_effect_fec = crop_effect_fec, dg = dg)
  seedbank_next = seedbank_post_germ + new_seeds 
  
  #plotting stuff
  max_below_ground = max(max(seedbank_current), max(seedbank_next))
  par(mfrow = c(1, 2))
  #above ground
  plot(eval_object$above_ground, new_plants[eval_object$above_ground_index], col = 'black', type = 'l', main = 'above ground', bty = 'n')
  lines(eval_object$above_ground, survivors_herb, col = 'red')
  lines(eval_object$above_ground, survivors_noherb, col = 'blue')
  lines(eval_object$above_ground, survivors_joint, col = 'green')
  legend(x = 'topright', legend = c('new_plants', 'survivors_herb', 'survivors_noherb', 'survival_joint'), col = c('black', 'red', 'blue', 'green'), lwd = 2)
  #put some lines at the mean to help see if dists have moved in expected direction
  mean_new_plants = get_mean_sd(new_plants[eval_object$above_ground_index], eval_object$above_ground, dg)$approx_mean
  lines(x = c(mean_new_plants, mean_new_plants), y = c(0, max(new_plants)), col = 'black')
  mean_survivors_herb = get_mean_sd(survivors_herb, eval_object$above_ground, dg)$approx_mean
  lines(x = c(mean_survivors_herb, mean_survivors_herb), y = c(0, max(survivors_herb)), col = 'red')
  mean_survivors_noherb = get_mean_sd(survivors_noherb, eval_object$above_ground, dg)$approx_mean
  lines(x = c(mean_survivors_noherb, mean_survivors_noherb), y = c(0, max(survivors_noherb)), col = 'blue')
  mean_survivors_joint = get_mean_sd(survivors_joint, eval_object$above_ground, dg)$approx_mean
  lines(x = c(mean_survivors_joint, mean_survivors_joint), y = c(0, max(survivors_joint)), col = 'green')
    
  #below ground
  plot(eval_object$seed, seq(0, max_below_ground, length = length(eval_object$seed)), type = 'n', bty = 'n', main = 'below ground')
  lines(eval_object$seed, seedbank_current, col = 'black')
  lines(eval_object$seed, new_seedbank, col = 'red')
  lines(eval_object$seed, seedbank_post_germ, col = 'blue')
  if(length(new_seeds) > 1){
    lines(eval_object$seed, new_seeds, col = 'green')
  }else{
    lines(eval_object$seed, rep(new_seeds, length(eval_object$seed)), col = 'green')
  }
  lines(eval_object$seed, seedbank_next, col = 'purple')
  legend(x = 'topright', legend = c('seedbank_current', 'new_seedbank', 'seedbank_post_germ', 'new_seeds', 'seedbank_next'), col = c('black', 'red', 'blue', 'green', 'purple'), lwd = 2)
  #put some lines at the mean to help see if dists have moved in expected direction
  mean_seedbank_current = get_mean_sd(seedbank_current, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_seedbank_current, mean_seedbank_current), y = c(0, max(seedbank_current)), col = 'black')
  mean_new_seedbank = get_mean_sd(new_seedbank, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_new_seedbank, mean_new_seedbank), y = c(0, max(new_seedbank)), col = 'red')
  mean_seedbank_post_germ = get_mean_sd(seedbank_post_germ, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_seedbank_post_germ, mean_seedbank_post_germ), y = c(0, max(seedbank_post_germ)), col = 'blue')
  mean_new_seeds = get_mean_sd(new_seeds, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_new_seeds, mean_new_seeds), y = c(0, max(new_seeds)), col = 'green')
  mean_seedbank_next = get_mean_sd(seedbank_next, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_seedbank_next, mean_seedbank_next), y = c(0, max(seedbank_next)), col = 'purple')
  
  #give the managment used on each turn 
  mtext(text = paste0('plowing: ', seed_movement, ' | herb: ', herb_rate, ' | mech: ', mech_control, ' | crop: ', crop_effect_sur, ' | density: ',  density_effect_fec), side = 3)
  
  #print out some numbers as we go along to get an idea of what is happening internally.
  print('###################################################################')
  print(paste0('###   iteration ', iteration, '   #####'))
  print(paste0('intial number seeds: ', sum(seedbank_current) * dg))
  print(paste0('num seeds after plow and seed sur: ', sum(new_seedbank) * dg))
  print(paste0('num seeds in soil post-germ: ', sum(seedbank_post_germ) * dg))
  print(paste0('num plants emerged post mech and crop sur: ', sum(new_plants) * dg))
  print(paste0('number herb survivors: ', sum(survivors_herb) * dg))
  print(paste0('number none herb surv: ', sum(survivors_noherb) * dg))
  print(paste0('total number sur (parents): ', sum(survivors_joint) * dg))
  print(paste0('num new seeds ', sum(new_seeds) * dg))
  print(paste0('num seeds next year ', sum(seedbank_next) * dg))
  
  return(seedbank_next)
}


##MULTI_ITERATION()
##plots multipule iterations for diagnostic purposes and to see what is happening 
multi_iteration_plot <- function(seedbank_initial, germination, seed_survival, eval_object, pro_exposed, sur0, sur_cost_resist, herb_effect, survive_resist, max_sur, fec_max, 
  fec0, fec_cost, offspring_sd, dense_depend_fec, dg, num_iter, sub_action, action_seq, output_loc, output_name = 'multi_iteration_output.pdf'){
  
  results = matrix(NA, nrow = num_iter, ncol = length(eval_object$seed))
  i = 1
  #create the function to do the offspring mixing 
  fecundity <- fecundity_closure(eval_points_object = eval_object, offspring_sd = offspring_sd, fec0 = fec0, fec_cost = fec_cost)
  intial_file_location = getwd()
  setwd(output_loc)
  pdf(file = output_name, width = 14, height = 7)
    results[i, ] = single_iteration_plot(seedbank_current = seedbank_initial, germination = germination, mech_control = sub_action$mech_sur[action_seq[i, 'mech']], 
    crop_effect_sur = sub_action$crop_sur[action_seq[i, 'crop']], seed_survival = seed_survival, seed_movement = sub_action$plow[action_seq[i, 'plow']], eval_object = eval_object, 
    pro_exposed = pro_exposed, herb_rate = sub_action$herb[action_seq[i, 'herb']], sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, 
    survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, dense_depend_fec = dense_depend_fec, crop_effect_fec = sub_action$crop_fec[action_seq[i, 'crop']], 
    density_effect_fec = sub_action$dens_fec[action_seq[i, 'dens']], fec_function = fecundity, iteration = i, dg = dg)

    for(i in 2:num_iter){
      results[i, ] = single_iteration_plot(seedbank_current = results[i - 1, ], germination = germination, mech_control = sub_action$mech_sur[action_seq[i, 'mech']], 
    crop_effect_sur = sub_action$crop_sur[action_seq[i, 'crop']], seed_survival = seed_survival, seed_movement = sub_action$plow[action_seq[i, 'plow']], eval_object = eval_object, 
    pro_exposed = pro_exposed, herb_rate = sub_action$herb[action_seq[i, 'herb']], sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, 
    survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, dense_depend_fec = dense_depend_fec, crop_effect_fec = sub_action$crop_fec[action_seq[i, 'crop']], 
    density_effect_fec = sub_action$dens_fec[action_seq[i, 'dens']], fec_function = fecundity, iteration = i, dg = dg)
    }
  dev.off()
  
  setwd(intial_file_location)
  
  results
}

#print version of the state update function
state_reward_next_econ_plot <- function(current_state, eval_object, action, sub_action, seed_survival, germination, pro_exposed, max_sur, sur0, sur_cost_resist, 
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
  new_seeds = fec_function(N_m = survivors_joint_mech, fec_max = fec_max, dense_depend_fec = dense_depend_fec, crop_effect_fec = sub_action$crop_fec[action['crop']], 
    density_effect_fec = sub_action$dens_fec[action['dens']], dg = dg)
  seedbank_next = seedbank_post_germ + new_seeds 
  #take these populations and new seedbank and turn them into rewards and states for the value function calculations 
  income = max(0, income0[action['crop']] - yeild_loss[action['crop']] * sum(survivors_joint_mech) * dg) #money made
  cost_nomech = cost_space_nomech$herb[action['herb']] + cost_space_nomech$crop[action['crop']] + cost_space_nomech$plow[action['plow']] + cost_space_nomech$dens[action['dens']] #non mech cost
  #mech cost
  if(action['mech'] == 1){
    cost_mech = 0
  }else{
    cost_mech = mech_cost0 + mech_cost * sum(survivors_joint) * dg
  }
  
  print('### parameters to survival calculation ###')
  print(c(max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist))
  print('eval_object$above_ground')
  print(eval_object$above_ground)
  print(c(herb = sub_action$herb[action['herb']], herb_effect = herb_effect, survive_resist = survive_resist))
  print(paste0('intial established pop: ', sum(new_plants) * dg))
  print(paste0('herbicide survivors: ', sum(survivors_herb) * dg))
  print(paste0('none herbicide survivors: ', sum(survivors_noherb) * dg))
  print(paste0('num_parents: ', sum(survivors_joint_mech) * dg))
  print(paste0('number new seeds: ', sum(new_seeds) * dg))
  print(paste0('num after none mech control: ', sum(survivors_joint) * dg))
  print(paste0('total yeild loss: ', yeild_loss[action['crop']] * sum(survivors_joint_mech) * dg))
  print(paste0('base income: ', income0[action['crop']]))
  print(paste0('net income: ', income))
  print(paste0('mech_cost0: ', mech_cost0))
  print(paste0('mech cost potential: ', mech_cost0 + mech_cost * sum(survivors_joint) * dg))
  print(paste0('mech cost spent: ', cost_mech))
  print(paste0('reward: ', income - cost_nomech - cost_mech))
  print('### parameter print END ###')
  
  #make a plot as a side effect
  domain_g = range(eval_object$above_ground)
  domain_pop = c(0, max(max(seedbank_current, seedbank_next)))
  plot(x = eval_object$seed, y = eval_object$seed, type = 'n', xlim = domain_g, ylim = domain_pop, bty = 'n', ylab = 'number ind')
  lines(eval_object$seed, seedbank_current, col = 'black')
  lines(eval_object$seed, new_seedbank, col = 'red')
  lines(eval_object$seed, seedbank_post_germ, col = 'blue')
  lines(eval_object$above_ground, survivors_herb, col = 'orange')
  lines(eval_object$above_ground, survivors_noherb, col = 'yellow')
  if(length(new_seeds) > 1){
    lines(eval_object$seed, new_seeds, col = 'green')
  }else{
    lines(eval_object$seed, rep(new_seeds, length(eval_object$seed)), col = 'green')
  }
  lines(eval_object$seed, seedbank_next, col = 'purple')
  legend(x = 'topright', legend = c('seedbank_current', 'new_seedbank', 'seedbank_post_germ', 'new_seeds', 'seedbank_next', 'survivors_herb', 'survivors_noherb'), 
    col = c('black', 'red', 'blue', 'green', 'purple', 'orange', 'yellow'), lwd = 2)
  #put some lines at the mean to help see if dists have moved in expected direction
  mean_seedbank_current = get_mean_sd(seedbank_current, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_seedbank_current, mean_seedbank_current), y = c(0, max(seedbank_current)), col = 'black')
  mean_new_seedbank = get_mean_sd(new_seedbank, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_new_seedbank, mean_new_seedbank), y = c(0, max(new_seedbank)), col = 'red')
  mean_seedbank_post_germ = get_mean_sd(seedbank_post_germ, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_seedbank_post_germ, mean_seedbank_post_germ), y = c(0, max(seedbank_post_germ)), col = 'blue')
  mean_new_seeds = get_mean_sd(new_seeds, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_new_seeds, mean_new_seeds), y = c(0, max(new_seeds)), col = 'green')
  mean_seedbank_next = get_mean_sd(seedbank_next, eval_object$seed, dg)$approx_mean
  lines(x = c(mean_seedbank_next, mean_seedbank_next), y = c(0, max(seedbank_next)), col = 'purple')
  mean_sur_herb = get_mean_sd(survivors_herb, eval_object$above_ground, dg)$approx_mean
  lines(x = c(mean_sur_herb, mean_sur_herb), y = c(0, max(mean_sur_herb)), col = 'orange')
  
  #give the managment used on each turn 
  mtext(text = paste0('plowing: ', sub_action$plow[action['plow']], ' | herb: ', sub_action$herb[action['herb']], ' | mech: ', sub_action$mech_sur[action['mech']],
    ' | crop_sur: ', sub_action$crop_sur[action['crop']], ' | density: ',  sub_action$dens_fec[action['dens']]), side = 3)
  
  return(list(current_state = current_state, next_state = build_state_1level(dist_top = seedbank_next, eval_points = eval_object$seed, dg = dg), reward = income - cost_nomech - cost_mech))
}

#simulate a policy from a set of starting points, record the state every time step, test various parts set it up so there is a bit more control over what actions are taken 
policy_simulator_test <- function(policy, model_pars, time_period, start_points, output_loc, output_name = 'policy_simulation_test.pdf'){
  mean_g_range = range(policy$sample_points_df$mean_g)
  pop_range = range(policy$sample_points_df$pop)
  #set up the scales training points for knn interpilation of the policy over state space
  mean_g_scaled = (policy$sample_points_df$mean_g - mean_g_range[1]) / (mean_g_range[2] - mean_g_range[1])
  pop_scaled = (policy$sample_points_df$pop - pop_range[1]) / (pop_range[2] - pop_range[1])
  #create fecundity function 
  fecundity = fecundity_closure(eval_points_object = model_pars$eval_obj, offspring_sd = model_pars$params['offspring_sd'], fec0 = model_pars$params['fec0'], 
    fec_cost = model_pars$params['fec_cost'])
    
  #print a vector of each parameter goint into state_reward_next_econ so can see if it is the expected one 
  #build the action vector
  new_point = data.frame(mean_g = (start_points[1, 'mean_g'] - mean_g_range[1]) / (mean_g_range[2] - mean_g_range[1]),
    pop = (start_points[1, 'pop'] - pop_range[1]) / (pop_range[2] - pop_range[1]))
  first_act = as.numeric(as.character(knn1(train = data.frame(mean_g = mean_g_scaled, pop = pop_scaled), test = new_point, cl = policy$sample_points_df$policy)))

  print('### parameters to simulate population ###')
  print(list(eval_object = model_pars$eval_obj, action = model_pars$action_space[first_act, ], sub_action = model_pars$sub_action, 
    seed_survival = model_pars$params['seed_survival'], germination = model_pars$params['germination'], pro_exposed = model_pars$params['pro_exposed'], 
    max_sur = model_pars$params['max_sur'], sur0 = model_pars$params['sur0'], sur_cost_resist = model_pars$params['sur_cost_resist'], herb_effect = model_pars$params['effect_herb'], 
    survive_resist = model_pars$params['survive_resist'], fec_max = model_pars$params['fec_max'], dense_depend_fec = model_pars$params['dense_depend_fec'], 
    income0 = model_pars$income0, yeild_loss = model_pars$yeild_loss, cost_space_nomech = model_pars$cost_space_nomech, mech_cost0 = model_pars$params['mech_cost0'], 
    mech_cost = model_pars$params['mech_cost'], fec_function = fecundity, dg = model_pars$params['dg'])) 
  print('##########################################')  
  
  #for every pointed sampled in the state space find the fate of the population at the time horizon for following that 
  int_loc = getwd()
  setwd(output_loc)
  pdf(file = output_name)
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
      #set the intial population USE THE PRINT VERSION OF THE 1 STEP SIMULATION TO SEE HOW THE POPS EVOLVE OVER TIME
      state_over_time[[1]] = state_reward_next_econ_plot(current_state = list(mean_top = s['mean_g'], sd_top = model_pars$params['pop_sd'], pop_top = s['pop']), 
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
	
	state_over_time[[ts]] = state_reward_next_econ_plot(current_state = state_over_time[[ts - 1]]$next_state, eval_object = model_pars$eval_obj, 
	  action = model_pars$action_space[action_seq[ts, 'act_num'], ], sub_action = model_pars$sub_action, seed_survival = model_pars$params['seed_survival'], 
	  germination = model_pars$params['germination'], pro_exposed = model_pars$params['pro_exposed'], max_sur = model_pars$params['max_sur'], sur0 = model_pars$params['sur0'], 
	  sur_cost_resist = model_pars$params['sur_cost_resist'], herb_effect = model_pars$params['effect_herb'], survive_resist = model_pars$params['survive_resist'], 
	  fec_max = model_pars$params['fec_max'], dense_depend_fec = model_pars$params['dense_depend_fec'], income0 = model_pars$income0, yeild_loss = model_pars$yeild_loss, 
	  cost_space_nomech = model_pars$cost_space_nomech, mech_cost0 = model_pars$params['mech_cost0'], mech_cost = model_pars$params['mech_cost'], fec_function = fecundity, 
	  dg = model_pars$params['dg'])
      }

      return(list(start_point = s, action_seq = action_seq, state_over_time = state_over_time))
    })  
    
  dev.off()  
  setwd(int_loc)
  return(results)
}



#Check the test results  
test_functions_broken <- function(file_loc){
  setwd(file_loc)
  test_obj_name <- load('nonspatial_model_test_answer_key.Rdata') #load the test key
  eval(parse(text = nonspace_test_answer_key[[1]]$question))#set parameters for the test run
  test1 = quant_gen_offspring_distribution(N_f = N_f, eval_points = eval_all$above_ground, offspring_sd = offspring_sd, seed_eval_points = eval_all$seed, dg = dg) #get the output form the current version of the function
  eval(parse(text = nonspace_test_answer_key[[2]]$question))#set parameters for the test run
  test2 = eval_points_builder(lower_eval_point = lower_eval_point, upper_eval_point = upper_eval_point, resolution = resolution, seed_expantion = seed_expantion)
  eval(parse(text = nonspace_test_answer_key[[3]]$question))#set parameters for the test run
  test3 = fecundity(N_m = N_m, eval_points = eval_points, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, N_f = N_f, offspring_sd = offspring_sd, seed_eval_points = seed_eval_points,
    dense_depend_fec = dense_depend_fec, crop_effect_fec = crop_effect_fec, density_effect_fec = density_effect_fec, dg = dg)
  eval(parse(text = nonspace_test_answer_key[[4]]$question))#set parameters for the test run
  test4 = single_iteration_1level(seedbank_current = seedbank_current, germination = germination, mech_control = mech_control, crop_effect_sur = crop_effect_sur, seed_survival = seed_survival, 
    seed_movement = seed_movement, eval_object = eval_object, pro_exposed = pro_exposed, herb_rate = herb_rate, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, 
    survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, 
    density_cutoff = density_cutoff, crop_effect_fec = crop_effect_fec, density_effect_fec = dense_depend_fec, dg = dg)
  eval(parse(text = nonspace_test_answer_key[[5]]$question))#set parameters for the test run
  test5 = eval_points_update(eval_points_object = eval_points_object, above_ground_dist = above_ground_dist, density_cutoff = density_cutoff)
  eval(parse(text = nonspace_test_answer_key[[6]]$question))#set parameters for the test run
  test6 = multi_iteration(seedbank_initial = seedbank_initial, germination = germination, seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, 
    sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost,
    offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, density_cutoff = density_cutoff, dg = dg, num_iter = num_iter, sub_action = sub_action, action_seq = action_seq)
  eval(parse(text = nonspace_test_answer_key[[7]]$question))
  test7 = state_reward_next_econ(current_state = current_state, eval_object = eval_object, action = action, sub_action = sub_action, seed_survival = seed_survival, 
    germination = germination, density_cutoff = density_cutoff, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, 
    survive_resist = survive_resist, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yield_loss, 
    cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, dg = dg)
  
  
  #test the seedbank functin shifts seeds without destroying them 
  x = seq(-50, 50, 0.1)
  test_sb0 = rbind(100 * dnorm(x, 10, 2.5), 20 * dnorm(x, 0, 2.5))
  test_sb = seedbank_plowing_2level(seedbank0 = test_sb0, seed_survival = 1, germination = 0, seed_movement = 0.2, deep_germ_reduce = 0)
  start_pop_top = sum(test_sb0[1, ] * 0.1)
  start_pop_bottom = sum(test_sb0[2, ] * 0.1)
  top_pop = sum(test_sb[1, ] * 0.1)
  bottom_pop = sum(test_sb[2, ] * 0.1)
  test_seedbank_num_constant = all.equal(top_pop + bottom_pop, start_pop_top + start_pop_bottom)
  test_seedbank_levels = all.equal(top_pop , 84) & all.equal(bottom_pop, 36)
  #test the state building functions and the state update function
  test_dist = 100 * dnorm(x, 10, 2.5)
  test_mean_sd_est = get_mean_sd(dist = test_dist, eval_points = x, dg = 0.1)
  test_equal_mean_sd = all.equal(test_mean_sd_est$approx_mean, 10) & all.equal(test_mean_sd_est$approx_sd, 2.5) & all.equal(test_mean_sd_est$total_pop, 100)
  
  
  test_results <- ifelse(all.equal(test1, nonspace_test_answer_key[[1]]$answer), 'QUANT_GEN_OFFSPRING_DISTRIBUTION() still fine', 'Something you did broke the function QUANT_GEN_OFFSPRING_DISTRIBUTION()')
  test_results[2] <- ifelse(identical(test2, nonspace_test_answer_key[[2]]$answer), 'EVAL_POINTS_BUILDER() still fine', 'Something you did broke the function EVAL_POINTS_BUILDER()')
  test_results[3] <- ifelse(identical(test3, nonspace_test_answer_key[[3]]$answer), 'FECUNDITY() still fine', 'Something you did broke the function FECUNDITY()')
  test_results[4] <- ifelse(identical(test4, nonspace_test_answer_key[[4]]$answer), 'SINGLE_INTERATION_1LEVEL() still fine', 'Something you did broke the function SINGLE_INTERATION_1LEVEL()') 
  test_results[5] <- ifelse(identical(test5, nonspace_test_answer_key[[5]]$answer), 'EVAL_POINTS_UPDATE() still fine', 'Something you did broke the function EVAL_POINTS_UPDATE()')
  test_results[6] <- ifelse(identical(test6, nonspace_test_answer_key[[6]]$answer), 'MULTI_ITERATION() still fine', 'Something you did broke the function MULTI_INTERATION()')
  test_results[7] <- ifelse(identical(test7, nonspace_test_answer_key[[7]]$answer), 'STATE_REWARD_NEXT_ECON() still fine', 'Something you did broke the function STATE_REWARD_NEXT_ECON()')
  test_results[8] <- ifelse(test_seedbank_num_constant & test_seedbank_levels, 'Two layer SEEDBANK() function numbers add up', 'Something broke the two layer SEEDBANK() and numbers no longer correct')
  test_results[9] <- ifelse(test_equal_mean_sd, 'GET_MEAN_SD() successfully recovers numbers', 'GET_MEAN_SD() cannot recover correct numbers')
  print(test_results) 
}


#add new questions and answers to the answer key 
#nonspace_test_answer_key[[7]] <- list()
#nonspace_test_answer_key[[7]]$question = paste0('eval_object = test2\ncurrent_state = build_state_1level(dist_top = dnorm(test2$seed, 5, 2.5), eval_points = test2$seed, dg = 0.5)\n',
#    'action = c(herb = 1, crop = 2, plow = 1, dens = 1, mech = 2)\nsub_action = list(herb = c(0, 1), crop_sur = c(1, 0.8, 0), mech_sur = c(1, 0.1), dens_fec = c(1, 0.9), crop_fec = c(1, 0.95, 1), plow = c(0, 0.3))\n',
#    'seed_survival = 0.5\ngermination = 0.8\ndensity_cutoff = 0.00001\npro_exposed = 0.8\nmax_sur = 0.95\nsur0 = 5\nsur_cost_resist = 0\nherb_effect = 5\nsurvive_resist = 0.8\n',
#    'fec_max = 100\nfec0 = 0\nfec_cost = 0.1\noffspring_sd = 0.7\ndense_depend_fec = 0.002\nincome0 = c(10, 2, 0)\nyield_loss = c(0.00001, 0.000001, 0)\n',
#    'cost_space_nomech = list(herb = c(0, 0.1), crop = c(0, 0.5, 0.5), plow = c(0, 0.05), dens = c(0, 0.5))\nmech_cost0 = 0.2\nmech_cost = 0.001\ndg = 0.5\n')
#eval(parse(text = nonspace_test_answer_key[[7]]$question))
#nonspace_test_answer_key[[7]]$answer = state_reward_next_econ(current_state = current_state, eval_object = eval_object, action = action, sub_action = sub_action, seed_survival = seed_survival, 
#  germination = germination, density_cutoff = density_cutoff, pro_exposed = pro_exposed, max_sur = max_sur, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, 
#  survive_resist = survive_resist, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, income0 = income0, yeild_loss = yield_loss, 
#  cost_space_nomech = cost_space_nomech, mech_cost0 = mech_cost0, mech_cost = mech_cost, dg = dg)
  
#save(nonspace_test_answer_key, file = 'nonspatial_model_test_answer_key.Rdata')
  

