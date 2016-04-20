#set of functions for herbicide resitance model covering various processes like genotype production and seed bank dynamics
#working_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_model' 
#setwd(working_loc)
#test_obj_name <- load('nonspatial_model_test_answer_key.Rdata') #load the test key

## EVAL_POINTS_BUILDER(lower_eval_point, upper_eval_point, resolution, seed_expantion)
## produces a matrix of evaluation points. The first row is evaluation points for evaluating above ground individuals, second row (which contains the first row) is evalution points for seeds
## lower_eval_point = lower value of g to be evaluated for above ground plants
## upper_eval_point = upper value of g to be evaluated for above ground plants
## resolution = resolution to evaluate distribtuions at
## seed_expantion = factor to spread the distribtuion by, typically will be some multipule of the sd of conditional breeding value distribution 
eval_points_builder <- function(lower_eval_point, upper_eval_point, resolution, seed_expantion){
  above_ground_eval = seq(lower_eval_point, upper_eval_point, resolution)
  seed_lower = seq(above_ground_eval[1] * seed_expantion, above_ground_eval[1], resolution)
  seed_upper = seq(above_ground_eval[length(above_ground_eval)], above_ground_eval[length(above_ground_eval)] * seed_expantion, resolution)
  seed_eval = c(seed_lower[1:(length(seed_lower) - 1)], above_ground_eval, seed_upper[2:length(seed_upper)])
  list(above_ground = above_ground_eval, seed = seed_eval, above_ground_index = which(seed_eval %in% above_ground_eval))
}

## EVAL_POINYS_UPDATE()
## Takes and eval_points_builder object and a above ground distribtuion of plants: returns an eval_points_builder object with an updated the above_ground eval window based on above ground distribution
##eval_points_object = names list produced by EVAL_POINTS_BUILDER() or a previous call to EVAL_POINTS_UPDATE()
## above_ground_dist = a distribution of above ground plants evaluated on seed eval points, 
## density_cutoff = population density (on distrbution over g) above which seed evaluation points are retained in the above ground evaluation points
eval_points_update <- function(eval_points_object, above_ground_dist, density_cutoff){
  eval_points_object$above_ground = eval_points_object$seed[which(above_ground_dist > density_cutoff)]
  if(length(eval_points_object$above_ground) < 10) eval_points_object$above_ground = eval_points_object$seed[which(above_ground_dist %in% tail(sort(above_ground_dist), 10))]
  eval_points_object$above_ground_index = which(eval_points_object$seed %in% eval_points_object$above_ground)
  eval_points_object
}

## produces a function which when called produces the mixing kernel of offsping over g given the parenrtal distribution 
##  offspring_sd is the sd of breeding value distribtuion of offspring (which we assume is normal).
#makes the closure to produce the function for the mixing kernel so a bunch of heavy calculation is pre-calculated and stored
quant_gen_closure <- function(eval_points_object, offspring_sd){
  ii <- expand.grid(M = 1:length(eval_points_object$above_ground), P = 1:length(eval_points_object$above_ground)) # index for every combination of maternal and paternal eval point 
  ep <- expand.grid(O = eval_points_object$seed, M = eval_points_object$above_ground, P = eval_points_object$above_ground)
  ep <- transform(ep, MP = 0.5 * M + 0.5 * P)[,c("O", "MP")]
  mix_kern <- dnorm(ep$O, mean = ep$MP, sd = offspring_sd)
  dim(mix_kern) <- c(length(eval_points_object$seed), length(eval_points_object$above_ground) ^ 2) #colSums(mix_kern) * dg = [1, 1, 1, ..., 1]
  rm(ep) 
 
  function(N_m){
    mix_kern %*% (N_m[ii$M] * N_m[ii$P])
  }
}

##FECUNDITY_CLOSURE()
## produces a function that when called returns the distribution of seeds on eval_points of g produced by population N_m (distribution of mothers on g evaluated at eval_points)
## N_m = maternal distrbution of indviduals over g
## eval_points = the values of g on which above ground individuals are evaluated
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## offspring_sd (passed to quant_gen_offspring_distribution()) = sd of conditional offspring distribution
## crop_effect_fec = proprional reduction in fecundity under different crops 
## density_effect_fec = reduction in density due to increased planting density 
## dense_depend_fec = 1/number of plants at which indivduals start to interfer with each other.
## dg = integration constant, equal to  eval_points[2] - eval_points[1]
#be careful this functions relies on N_m being evaluated on eval_points before being passed to this function so that the indexes all match up
fecundity_closure <- function(eval_points_object, offspring_sd, fec0, fec_cost){
  #constants for the mixing distribution
  ii <- expand.grid(M = 1:length(eval_points_object$above_ground), P = 1:length(eval_points_object$above_ground)) # index for every combination of maternal and paternal eval point 
  ep <- expand.grid(O = eval_points_object$seed, M = eval_points_object$above_ground, P = eval_points_object$above_ground)
  ep <- transform(ep, MP = 0.5 * M + 0.5 * P)[,c("O", "MP")]
  mix_kern <- dnorm(ep$O, mean = ep$MP, sd = offspring_sd)
  dim(mix_kern) <- c(length(eval_points_object$seed), length(eval_points_object$above_ground) ^ 2) #colSums(mix_kern) * dg = [1, 1, 1, ..., 1]
  rm(ep) 
  #constants for the fecundity function 
  resist_effect_fec = exp(-(fec0 - eval_points_object$above_ground * fec_cost))
  
  function(N_m, fec_max, dense_depend_fec, dg){
    num_survivors = sum(N_m) * dg
    if(num_survivors > 0){
      repo_happiness = 1 / (1 + resist_effect_fec + dense_depend_fec * num_survivors + dense_depend_fec * num_survivors * resist_effect_fec)
      post_select_parents = N_m * repo_happiness
      post_select_parents_normz = post_select_parents / (sum(post_select_parents) * dg)
      pd_seeds = (mix_kern %*% (post_select_parents_normz[ii$M] * post_select_parents_normz[ii$P])) * dg * dg #need to integrate 2 times across mix_kern cols so whole thing intergates to 1
      return(sum(repo_happiness * post_select_parents) * fec_max * dg * pd_seeds) #some slight inaccuracy here, if the expected number of total seeds is 20,000 the actual total number is 20,0047 when dg = 0.5, can be made much more accurate if dg is smaller  
    }else{
      return(0)
    }
  }
}


## SEEDBANK_PLOWING_2LEVEL()
## produces a distribution of seeds in the seed bank over the eval_points on g. 
## seedbank0 = distrbution of seeds in the seedbank in the last timestep
## seed_survival = probability that a seed in the seed bank survives one year
## germination = the probability that a seed in the seedbank survies one timestep
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## new_seeds = distribution of seeds produced by call to FECUNDITY()
## seed_movement = the proportion of seeds that get moved from one level to another
## deep_germ_reduce = the proportional reduction in germination for seeds in the bottom level of the seedbank. 
seedbank_plowing_2level <- function(seedbank0, seed_survival, germination, seed_movement, deep_germ_reduce){
  seedbank_top = (seedbank0[1, ] - seedbank0[1, ] *  seed_movement + seedbank0[2, ] * seed_movement) * seed_survival * (1 - germination)
  seedbank_bottom = (seedbank0[2, ] - seedbank0[2, ] *  seed_movement + seedbank0[1, ] * seed_movement)* seed_survival * (1 - germination * deep_germ_reduce)
  res = rbind(seedbank_top, seedbank_bottom)
  return(res)
}

## NEW_SEEDS()
## adds seeds produced by FECUNDITY() to the seedbank
## seedbank = two layer seed bank produced by call to SEEDBANK_PLOWING()
## newseeds = distribution of seeds produced by call to FECUNDITY()
new_seeds <- function(seedbank, newseeds){
  updated_seedbank = seedbank
  updated_seedbank[1, ] = updated_seedbank[1, ] + newseeds
  return(updated_seedbank)
}

## SINGLE_INTERATION_1LEVEL()
## produces a distrbution over g of indivduals in the seed bank including survival, reproduction and emergence 
## seedbank_current = distribution of seeds in the seed bank over g (either an intial seedbank or returned from another call to SINGLE_INTERATION()) 
## germination = germination probability
## pro_exposed = proprtion of emerged individuals that are exposed to herbicide
## seed_survival = probability that a seed in the seed bank survives one year
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## fec_max = the maximum number of seeds per mothers
## dense_depend_fec = density dependent effect on fecundity in units of 1/number individuals at which they start to affect each others fecundity
## herb_rate = 0 or 1 factor that if herbicide was applied or not
## sur0 = susrvial rate when g is 0 and there is no herbicide
## herb_effect = effect of herbicide on survival
## survive_resist = protective effect of a one unit increase in resistance score g
## dg = difference between evaluation points, should = eval_points[2] - eval_points[1]

single_iteration_1level <- function(seedbank_current, germination, seed_survival, eval_object, pro_exposed, herb_rate, survive_resist, 
  sur0, sur_cost_resist, herb_effect, fec_max, dense_depend_fec, fec_function, dg){

  new_seedbank = seedbank_current * seed_survival 
  new_plants = new_seedbank * germination 
  seedbank_post_germ = new_seedbank * (1 - germination)
  survivors_herb = pro_exposed * new_plants[eval_object$above_ground_index] *  (1 / (1 + exp(-(sur0 - herb_rate * (herb_effect - 
    pmin(herb_effect, survive_resist * eval_object$above_ground))))))
  survivors_noherb = (1 - pro_exposed) * new_plants[eval_object$above_ground_index] *  (1 / (1 + exp(-(sur0))))
  survivors_joint = survivors_herb + survivors_noherb
  seedbank_next = seedbank_post_germ + fec_function(N_m = survivors_joint, fec_max = fec_max, dense_depend_fec = dense_depend_fec, 
    dg = dg) 
  
  return(list(seedbank = seedbank_next, above_ground_pre_herb = new_plants[eval_object$above_ground_index], above_ground_post_herb = survivors_joint))
}

##UPDATE_TSR()
## takes a named vector of TSR genotypes (RR, Rr, rr) where each element is the current propotion of seedbank in each genotype
## a vector of strings giving the resistiant genotypes (must be in names(TSR_pop_current))
## s0, base survival at g = 0 no herb
## pro_exposed to herbiced
## number of above ground platns post herbicide
## number of above ground plants pre-herbicide 
## number of seeds in the seedbank currently (post germination)
## number of new seeds produced in the year
update_TSR <- function(TSR_freq_current, res_genotypes, sur0, pro_exposed, num_post_herb, 
  num_pre_herb, num_seedbank, num_new_seed){

  # proportions after survival
  TSR_post_sur = TSR_freq_current * ifelse(names(TSR_freq_current) %in% res_genotypes, 
  1 / (1 + exp(-sur0)), #TS resistant survival
  (1 - pro_exposed) / (1 + exp(-sur0)) + pro_exposed * (num_post_herb / num_pre_herb)) #TS non-resistant survival
  
  # frequency of each allele in the surviving population 
  freq_R = as.numeric((2 * TSR_post_sur['RR'] + TSR_post_sur['Rr']) / (2 * sum(TSR_post_sur)))
  freq_r = as.numeric((2 * TSR_post_sur['rr'] + TSR_post_sur['Rr']) / (2 * sum(TSR_post_sur)))
  
  freq_G_new_seed = c(RR = freq_R ^ 2, Rr = 2 * freq_R * freq_r, rr = freq_r ^2)
 
  total_seed = TSR_freq_current * num_seedbank + freq_G_new_seed * num_new_seed
  old_v_new = ifelse(total_seed > 0, (TSR_freq_current * num_seedbank) / total_seed, 0)
    
  TSR_freq_next = (old_v_new * TSR_freq_current + (1 - old_v_new) * freq_G_new_seed) / 
    (sum(old_v_new * TSR_freq_current + (1 - old_v_new) * freq_G_new_seed))

  return(TSR_freq_next)
}

## SINGLE_ITERATION_TSR()
## produces a distirbution of indivduals over g when there is both metabolic and TS resistance
## TSR_pop_current a named vector where each element is the proportion of population for one of thee res_genotypes
## res_genotypes, the genotypes in names(TSR_pop_current) that give target site resistance.
single_iteration_TSR <- function(seedbank_current, germination, seed_survival, eval_object, pro_exposed, herb_rate, survive_resist, 
  sur0, herb_effect, fec_max, dense_depend_fec, TSR_freq_current, res_genotypes, fec_function, dg){

  new_seedbank = seedbank_current * seed_survival 
  new_plants = new_seedbank * germination 
  seedbank_post_germ = new_seedbank * (1 - germination)
  pro_not_affected = (1 - pro_exposed) + pro_exposed * TSR_freq_current[res_genotypes]
  survivors_herb = new_plants[eval_object$above_ground_index] / (1 + exp(-(sur0 - herb_rate * (herb_effect - 
    pmin(herb_effect, survive_resist * eval_object$above_ground)))))
  survivors_noherb = new_plants[eval_object$above_ground_index]  / (1 + exp(-(sur0)))
  survivors_joint = (1 - pro_not_affected) * survivors_herb +  pro_not_affected * survivors_noherb
 
  new_seeds = fec_function(N_m = survivors_joint, fec_max = fec_max, dense_depend_fec = dense_depend_fec, dg = dg)
  seedbank_next = seedbank_post_germ + new_seeds 
    
  #track level of each TSR genotype
  num_post_herb = sum(survivors_herb) * dg 
  num_pre_herb = sum(new_plants) * dg
  num_seedbank = sum(seedbank_post_germ) * dg 
  num_new_seed = sum(new_seeds) * dg
    
  TSR_freq_next = update_TSR(TSR_freq_current = TSR_freq_current, res_genotypes = res_genotypes, 
    sur0 = sur0, pro_exposed = pro_exposed, num_post_herb = num_post_herb, num_pre_herb = num_pre_herb, 
    num_seedbank = num_seedbank, num_new_seed = num_new_seed)
    
  return(list(seedbank = seedbank_next, above_ground_pre_herb = new_plants[eval_object$above_ground_index], 
    above_ground_post_herb = survivors_joint, TSR_freq_next = TSR_freq_next))
}

## MULTI_ITERATION_TSR()
## produces a num_iter by length(eval_object$seed) matrix where each row is a distrbution over g of indivduals in the seedbank for 1 iteration   
## seedbank_current = distribution of seeds in the seed bank over g (either an intial seedbank or returned from another call to SINGLE_INTERATION()) 
## germination = germination probability
## pro_exposed = proprtion of emerged individuals that are exposed to herbicide
## seed_survival = probability that a seed in the seed bank survives one year
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## offspring_sd (passed to quant_gen_offspring_distribution()) = sd of conditional offspring distribution
## dense_depend_fec = density dependent effect on fecundity in units of 1/number individuals at which they start to affect each others fecundity
## herb_rate = 0 or 1 factor that if herbicide was applied or not
## sur0 = susrvial rate when g is 0 and there is no herbicide
## herb_effect = effect of herbicide on survival
## survive_resist = protective effect of a one unit increase in resistance score g
## dg = difference between evaluation points, should = eval_points[2] - eval_points[1]
## num_iter = number of interations to run the model for
## action_seq = matrix 5 x num_iter matrix that gives the sub-action taken in every iteration, used to reference action_space
multi_iteration_TSR <- function(seedbank_initial, TSR_inital, germination, seed_survival, 
  eval_object, pro_exposed, sur0, herb_effect, survive_resist, fec_max, fec0, fec_cost, 
  offspring_sd, dense_depend_fec, dg, num_iter, herb_seq, res_genotypes){
    
  seedbank_results = matrix(NA, nrow = num_iter, ncol = length(eval_object$seed))
  ag_pre_herb_results = matrix(NA, nrow = num_iter, ncol = length(eval_object$above_ground))
  ag_post_herb_results = matrix(NA, nrow = num_iter, ncol = length(eval_object$above_ground))
  TSR_results = matrix(NA, nrow = num_iter, ncol  = length(TSR_inital))
  colnames(TSR_results) <- names(TSR_inital)
  
  #create the function to do the offspring mixing 
  fecundity <- fecundity_closure(eval_points_object = eval_object, offspring_sd = offspring_sd, 
    fec0 = fec0, fec_cost = fec_cost)
  i = 1
  results = single_iteration_TSR(seedbank_current = seedbank_initial, germination = germination,  
    seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, herb_rate = herb_seq[i], 
    TSR_freq_current = TSR_inital, res_genotypes = res_genotypes, sur0 = sur0, herb_effect = herb_effect, 
    survive_resist = survive_resist, fec_max = fec_max, dense_depend_fec = dense_depend_fec, 
    fec_function = fecundity, dg = dg)
    
  seedbank_results[i, ] = results$seedbank
  ag_pre_herb_results[i, ] = results$above_ground_pre_herb
  ag_post_herb_results[i, ] = results$above_ground_post_herb
  TSR_results[i, ] = results$TSR_freq_next
  
  for(i in 2:num_iter){
    results = single_iteration_TSR(seedbank_current = seedbank_results[i - 1, ], germination = germination,  
      seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, 
      herb_rate = herb_seq[i], TSR_freq_current = TSR_results[i - 1, ], res_genotypes = res_genotypes, 
      sur0 = sur0, herb_effect = herb_effect, survive_resist = survive_resist, fec_max = fec_max, 
      dense_depend_fec = dense_depend_fec, fec_function = fecundity, dg = dg)
      
    seedbank_results[i, ] = results$seedbank
    ag_pre_herb_results[i, ] = results$above_ground_pre_herb
    ag_post_herb_results[i, ] = results$above_ground_post_herb
    TSR_results[i, ] = results$TSR_freq_next
  }
  
  return(list(seedbank = seedbank_results, above_ground_pre_herb = ag_pre_herb_results, 
    above_ground_post_herb = ag_post_herb_results, TSR_freq = TSR_results))
}

## POP_RUN_TSR()
## runs a set of parameters under both herb and no_herb, after a burn in period
pop_run_TSR <- function(param_vect, int_pop_size, int_sd_g, res_genotypes, sur0 = 10,
  offspring_sd = 1, burnin_time, num_iter, eval_object, dg){
  
  #build an intial population
  int_pop = dnorm(eval_object$seed, 0, int_sd_g) * int_pop_size
  
  #run this population for a burnin period to get sd and mean g of a niave population given the parameter values  
  burnin_pop = multi_iteration_TSR(seedbank_initial = int_pop, TSR_inital = c(RR = 0, Rr = 0, rr = 1), germination = param_vect['germ_prob'], 
    seed_survival = param_vect['seed_sur'], eval_object = eval_object, pro_exposed = param_vect['pro_exposed'], 
    sur0 = sur0, herb_effect = param_vect['herb_effect_mult'] * sur0, survive_resist = sur0 * param_vect['herb_effect_mult'] * param_vect['sur_protect_mult'], 
    fec_max = param_vect['fec_max'], fec0 = param_vect['fec0'], fec_cost = param_vect['fec_cost_mult'] * param_vect['fec0'], 
    offspring_sd = offspring_sd, dense_depend_fec = param_vect['fec_dd'], dg = dg, num_iter = burnin_time, 
    herb_seq = rep(0, burnin_time), res_genotypes) #build the population 
  
  niave_mean_sd = get_mean_sd(dist = burnin_pop$seedbank[burnin_time, ], eval_points = eval_object$seed, dg = dg)
  niave_pop = dnorm(eval_object$seed, niave_mean_sd$approx_mean, niave_mean_sd$approx_sd) * int_pop_size
  
  intial_Rr = as.numeric(param_vect['int_Rr'])
  #run the population both with and without herbicide 
  noherb_pop = multi_iteration_TSR(seedbank_initial = niave_pop, TSR_inital = c(RR = 0, Rr = intial_Rr, rr = 1 - intial_Rr), 
    germination = param_vect['germ_prob'], seed_survival = param_vect['seed_sur'], eval_object = eval_object, pro_exposed = param_vect['pro_exposed'], 
    sur0 = sur0, herb_effect = param_vect['herb_effect_mult'] * sur0, survive_resist = sur0 * param_vect['herb_effect_mult'] * param_vect['sur_protect_mult'], 
    fec_max = param_vect['fec_max'], fec0 = param_vect['fec0'], fec_cost = param_vect['fec_cost_mult'] * param_vect['fec0'], 
    offspring_sd = offspring_sd, dense_depend_fec = param_vect['fec_dd'], dg = dg, num_iter = num_iter, 
    herb_seq = rep(0, num_iter), res_genotypes)
    
  herb_pop = multi_iteration_TSR(seedbank_initial = niave_pop, TSR_inital = c(RR = 0, Rr = intial_Rr, rr = 1 - intial_Rr), 
    germination = param_vect['germ_prob'], seed_survival = param_vect['seed_sur'], eval_object = eval_object, pro_exposed = param_vect['pro_exposed'], 
    sur0 = sur0, herb_effect = param_vect['herb_effect_mult'] * sur0, survive_resist = sur0 * param_vect['herb_effect_mult'] * param_vect['sur_protect_mult'], 
    fec_max = param_vect['fec_max'], fec0 = param_vect['fec0'], fec_cost = param_vect['fec_cost_mult'] * param_vect['fec0'], 
    offspring_sd = offspring_sd, dense_depend_fec = param_vect['fec_dd'], dg = dg, num_iter = num_iter, 
    herb_seq = rep(1, num_iter), res_genotypes)

  return(list(parameters = param_vect, burnin_pop = burnin_pop, noherb_pop = noherb_pop, herb_pop = herb_pop)) 
}


## MULTI_ITERATION()
## produces a num_iter by length(eval_object$seed) matrix where each row is a distrbution over g of indivduals in the seedbank for 1 iteration   
## seedbank_current = distribution of seeds in the seed bank over g (either an intial seedbank or returned from another call to SINGLE_INTERATION()) 
## germination = germination probability
## pro_exposed = proprtion of emerged individuals that are exposed to herbicide
## seed_survival = probability that a seed in the seed bank survives one year
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## offspring_sd (passed to quant_gen_offspring_distribution()) = sd of conditional offspring distribution
## dense_depend_fec = density dependent effect on fecundity in units of 1/number individuals at which they start to affect each others fecundity
## herb_rate = 0 or 1 factor that if herbicide was applied or not
## sur0 = susrvial rate when g is 0 and there is no herbicide
## herb_effect = effect of herbicide on survival
## survive_resist = protective effect of a one unit increase in resistance score g
## dg = difference between evaluation points, should = eval_points[2] - eval_points[1]
## num_iter = number of interations to run the model for
## action_seq = matrix 5 x num_iter matrix that gives the sub-action taken in every iteration, used to reference action_space
multi_iteration <- function(seedbank_initial, germination, seed_survival, eval_object, pro_exposed, sur0, herb_effect, survive_resist, 
  fec_max, fec0, fec_cost, offspring_sd, dense_depend_fec, dg, num_iter, herb_seq){
    
  seedbank_results = matrix(NA, nrow = num_iter, ncol = length(eval_object$seed))
  ag_pre_herb_results = matrix(NA, nrow = num_iter, ncol = length(eval_object$above_ground))
  ag_post_herb_results = matrix(NA, nrow = num_iter, ncol = length(eval_object$above_ground))
  
  #create the function to do the offspring mixing 
  fecundity <- fecundity_closure(eval_points_object = eval_object, offspring_sd = offspring_sd, fec0 = fec0, fec_cost = fec_cost)
  i = 1
  results = single_iteration_1level(seedbank_current = seedbank_initial, germination = germination,  
     seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, herb_rate = herb_seq[i], 
     sur0 = sur0, herb_effect = herb_effect, survive_resist = survive_resist, fec_max = fec_max, 
     dense_depend_fec = dense_depend_fec, fec_function = fecundity, dg = dg)
  seedbank_results[i, ] = results$seedbank
  ag_pre_herb_results[i, ] = results$above_ground_pre_herb
  ag_post_herb_results[i, ] = results$above_ground_post_herb
  
  for(i in 2:num_iter){
    results = single_iteration_1level(seedbank_current = seedbank_results[i - 1, ], germination = germination,  
      seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, herb_rate = herb_seq[i], 
      sur0 = sur0, herb_effect = herb_effect, survive_resist = survive_resist, fec_max = fec_max, dense_depend_fec = dense_depend_fec, 
      fec_function = fecundity, dg = dg)
      
    seedbank_results[i, ] = results$seedbank
    ag_pre_herb_results[i, ] = results$above_ground_pre_herb
    ag_post_herb_results[i, ] = results$above_ground_post_herb
  }
  
  return(list(seedbank = seedbank_results, above_ground_pre_herb = ag_pre_herb_results, above_ground_post_herb = ag_post_herb_results))
}

## GET_EQUIL_SD()
## funtion to find the equilibrium sd for a given population size.
## found by iterating the population, but resetting the population size to a fixed
## at each time step
get_equil_sd <- function(pop_size, mean_g, int_sd, germination, seed_survival, eval_object, 
  pro_exposed, sur0, herb_effect, survive_resist, fec_max, fec0, fec_cost, offspring_sd, 
  dense_depend_fec, dg, num_iter){
  # set up list to hold the results  
  result = list() 
  #create the function to do the offspring mixing 
  fecundity <- fecundity_closure(eval_points_object = eval_object, offspring_sd = offspring_sd, fec0 = fec0, fec_cost = fec_cost)
  
  #make every iteration no herb
  herb_seq = rep(0, num_iter)
  
  #buiild inital population
  old_pop = dnorm(eval_object$seed, mean_g, int_sd) * pop_size
  
  #after each iteration take the mean and sd of the resulting population and build new pop of desired size
  for(i in 1:num_iter){
    new_pop = single_iteration_1level(seedbank_current = old_pop, germination = germination,  
      seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, herb_rate = herb_seq[i], 
      sur0 = sur0, herb_effect = herb_effect, survive_resist = survive_resist, fec_max = fec_max, 
      dense_depend_fec = dense_depend_fec, fec_function = fecundity, dg = dg)$seedbank
     
    result[[i]] = get_mean_sd(dist = new_pop, eval_points = eval_object$seed, dg = dg)
     
    #build the old pop for the next go round
    old_pop = dnorm(eval_object$seed, mean_g, result[[i]]$approx_sd) * pop_size
  
  }
 
  return(list(sd_g = sapply(result, FUN = function(x) x$approx_sd), pop = sapply(result, FUN = function(x) x$total_pop),
    mean_g = sapply(result, FUN = function(x) x$approx_mean)))
}

## function to pass to mapply to do the parameter sweeps. 
pop_run <- function(param_vect, int_pop_size, int_sd_g, sur0 = 10, offspring_sd = 1, burnin_time, num_iter, eval_object, dg){
  
   #build an intial population
  int_pop = dnorm(eval_object$seed, 0, int_sd_g) * int_pop_size
  
  #run this population for a burnin period to get sd and mean g of a niave population given the parameter values  
  burnin_pop = multi_iteration(seedbank_initial = int_pop, germination = param_vect['germ_prob'], 
    seed_survival = param_vect['seed_sur'], eval_object = eval_object, pro_exposed = param_vect['pro_exposed'], 
    sur0 = sur0, herb_effect = param_vect['herb_effect_mult'] * sur0, survive_resist = sur0 * param_vect['herb_effect_mult'] * param_vect['sur_protect_mult'], 
    fec_max = param_vect['fec_max'], fec0 = param_vect['fec0'], fec_cost = param_vect['fec_cost_mult'] * param_vect['fec0'], 
    offspring_sd = offspring_sd, dense_depend_fec = param_vect['fec_dd'], dg = dg, num_iter = burnin_time, 
    herb_seq = rep(0, burnin_time)) #build the population 
  
  niave_mean_sd = get_mean_sd(dist = burnin_pop$seedbank[burnin_time, ], eval_points = eval_object$seed, dg = dg)
  niave_pop = dnorm(eval_object$seed, niave_mean_sd$approx_mean, niave_mean_sd$approx_sd) * int_pop_size
  #run the population both with and without herbicide 
  noherb_pop = multi_iteration(seedbank_initial = niave_pop, germination = param_vect['germ_prob'],
    seed_survival = param_vect['seed_sur'], eval_object = eval_object, pro_exposed = param_vect['pro_exposed'], 
    sur0 = sur0, herb_effect = param_vect['herb_effect_mult'] * sur0, survive_resist = sur0 * param_vect['herb_effect_mult'] * param_vect['sur_protect_mult'], 
    fec_max = param_vect['fec_max'], fec0 = param_vect['fec0'], fec_cost = param_vect['fec_cost_mult'] * param_vect['fec0'], 
    offspring_sd = offspring_sd, dense_depend_fec = param_vect['fec_dd'], dg = dg, num_iter = num_iter, 
    herb_seq = rep(0, num_iter))
  herb_pop = multi_iteration(seedbank_initial = niave_pop, germination = param_vect['germ_prob'], 
    seed_survival = param_vect['seed_sur'], eval_object = eval_object, pro_exposed = param_vect['pro_exposed'], 
    sur0 = sur0, herb_effect = param_vect['herb_effect_mult'] * sur0, survive_resist = sur0 * param_vect['herb_effect_mult'] * param_vect['sur_protect_mult'], 
    fec_max = param_vect['fec_max'], fec0 = param_vect['fec0'], fec_cost = param_vect['fec_cost_mult'] * param_vect['fec0'], 
    offspring_sd = offspring_sd, dense_depend_fec = param_vect['fec_dd'], dg = dg, num_iter = num_iter, 
    herb_seq = rep(1, num_iter))

  return(list(parameters = param_vect, burnin_pop = burnin_pop, noherb_pop = noherb_pop, herb_pop = herb_pop)) 
}


## function that takes a list of population matricies, one with herb, one without, and a vector of parameters,
## returns a vector of population sizes over time with and without herb along with parameter vector
pops_list_cleanup <- function(pops_list){
  
  param_mat = sapply(pops_list, FUN = function(x) x$parameters)
 
  TSR_freq_list = lapply(pops_list, FUN = function(x) x$herb_pop$TSR_freq)
  g_seedbank_herb = sapply(pops_list, FUN = function(x){
    sapply(apply(x$herb_pop$seedbank, MARGIN = 1, FUN = get_mean_sd, eval_points = eval_object$seed, dg = dg), 
      FUN = function(x) x$approx_mean)
  })
  pop_seedbank_herb = sapply(pops_list, FUN = function(x){
    sapply(apply(x$herb_pop$seedbank, MARGIN = 1, FUN = get_mean_sd, eval_points = eval_object$seed, dg = dg), 
      FUN = function(x) x$total_pop)
  })
  pop_ag_post_herb = sapply(pops_list, FUN = function(x){
    sapply(apply(x$herb_pop$above_ground_post_herb, MARGIN = 1, FUN = get_mean_sd, eval_points = eval_object$above_ground, dg = dg), 
      FUN = function(x) x$total_pop)
  })
  pop_ag_pre_herb = sapply(pops_list, FUN = function(x){
    sapply(apply(x$herb_pop$above_ground_pre_herb, MARGIN = 1, FUN = get_mean_sd, eval_points = eval_object$above_ground, dg = dg), 
      FUN = function(x) x$total_pop)
  })
  pop_seedbank_noherb = sapply(pops_list, FUN = function(x){
    sapply(apply(x$noherb_pop$seedbank, MARGIN = 1, FUN = get_mean_sd, eval_points = eval_object$seed, dg = dg), 
      FUN = function(x) x$total_pop)
  })
  pop_ag_noherb = sapply(pops_list, FUN = function(x){
    sapply(apply(x$noherb_pop$above_ground_post_herb, MARGIN = 1, FUN = get_mean_sd, eval_points = eval_object$above_ground, dg = dg), 
      FUN = function(x) x$total_pop)
  })
  
  return(list(parameters = param_mat, seedbank_herb = pop_seedbank_herb, seedbank_noherb = pop_seedbank_noherb,
    ag_pre_herb = pop_ag_pre_herb, ag_post_herb = pop_ag_post_herb, ag_noherb = pop_ag_noherb, TSR_freq = TSR_freq_list, 
    seedbank_g_herb = g_seedbank_herb))
}


## takes the output of pop_list_cleanup and calulates a set of metrics of population performance
## loking at max populations, final populations and comparisons of herb vs no_herb pops
pop_metrics <- function(pop_mats){

  #seedbank
  sb_over_yield = sapply(1:dim(pop_mats$seedbank_herb)[2], 
    FUN = function(x) sum(pop_mats$seedbank_herb[, x] > pop_mats$seedbank_noherb[, x]) >= 1)
  time_2_oy_sb = sapply(1:dim(pop_mats$seedbank_herb)[2], 
    FUN = function(x){ 
      if(sb_over_yield[x]){
	return(which(pop_mats$seedbank_herb[, x] > pop_mats$seedbank_noherb[, x])[1])
      }else{
	return(NA)
      }
    }) 
  time_in_oy_sb = sapply(1:dim(pop_mats$seedbank_herb)[2], 
    FUN = function(x){ 
      if(sb_over_yield[x]){
	return(sum(pop_mats$seedbank_herb[, x] > pop_mats$seedbank_noherb[, x]))
      }else{
	return(NA)
      }
    }) 
  mag_oy_sb = sapply(1:dim(pop_mats$seedbank_herb)[2], 
    FUN = function(x){ 
      if(sb_over_yield[x]){
	oy_elements_herb = pop_mats$seedbank_herb[pop_mats$seedbank_herb[, x] > pop_mats$seedbank_noherb[, x], ]
	oy_elements_noherb = pop_mats$seedbank_noherb[pop_mats$seedbank_herb[, x] > pop_mats$seedbank_noherb[, x], ]
	return(max(oy_elements_herb - oy_elements_noherb))
      }else{
	return(NA)
      }
    }) 
  diff_sb = sapply(1:dim(pop_mats$seedbank_herb)[2], 
    FUN = function(x){ 
      min(pop_mats$seedbank_noherb[, x] - pop_mats$seedbank_herb[, x])
    }) 
 
  #above ground
  ag_over_yield = sapply(1:dim(pop_mats$ag_herb)[2], 
    FUN = function(x) sum(pop_mats$ag_herb[, x] > pop_mats$ag_noherb[, x]) >= 1)
  time_2_oy_ag = sapply(1:dim(pop_mats$ag_herb)[2], 
    FUN = function(x){ 
      if(ag_over_yield[x]){
	return(which(pop_mats$ag_herb[, x] > pop_mats$ag_noherb[, x])[1])
      }else{
	return(NA)
      }
    }) 
  time_in_oy_ag = sapply(1:dim(pop_mats$ag_herb)[2], 
    FUN = function(x){ 
      if(ag_over_yield[x]){
	return(sum(pop_mats$ag_herb[, x] > pop_mats$ag_noherb[, x]))
      }else{
	return(NA)
      }
    }) 
  mag_oy_ag = sapply(1:dim(pop_mats$ag_herb)[2], 
    FUN = function(x){ 
      if(ag_over_yield[x]){
	oy_elements_herb = pop_mats$ag_herb[pop_mats$ag_herb[, x] > pop_mats$ag_noherb[, x], ]
	oy_elements_noherb = pop_mats$ag_noherb[pop_mats$ag_herb[, x] > pop_mats$ag_noherb[, x], ]
	return(max(oy_elements_herb - oy_elements_noherb))
      }else{
	return(NA)
      }
    }) 
  diff_ag = sapply(1:dim(pop_mats$seedbank_herb)[2], 
    FUN = function(x){ 
      min(pop_mats$ag_noherb[, x] - pop_mats$ag_herb[, x])
    }) 
   
  # some over all metrics of each line 
  fin_ag_noherb = pop_mats$ag_noherb[dim(pop_mats$ag_noherb)[1], ]
  fin_ag_herb = pop_mats$ag_herb[dim(pop_mats$ag_herb)[1], ]
  max_ag_noherb = apply(pop_mats$ag_noherb, MARGIN = 2, FUN = max)
  max_ag_herb = apply(pop_mats$ag_herb, MARGIN = 2, FUN = max)
  time_2_max_ag_noherb = apply(pop_mats$ag_noherb, MARGIN = 2, FUN = function(x) which(x == max(x))[1])
  time_2_max_ag_herb = apply(pop_mats$ag_herb, MARGIN = 2, FUN = function(x) which(x == max(x))[1])
  fin_sb_noherb = pop_mats$seedbank_noherb[dim(pop_mats$seedbank_noherb)[1], ]
  fin_sb_herb = pop_mats$seedbank_herb[dim(pop_mats$seedbank_herb)[1], ]
  max_sb_noherb = apply(pop_mats$seedbank_noherb, MARGIN = 2, FUN = max)
  max_sb_herb = apply(pop_mats$seedbank_herb, MARGIN = 2, FUN = max)
  time_2_max_sb_noherb = apply(pop_mats$seedbank_noherb, MARGIN = 2, FUN = function(x) which(x == max(x))[1])
  time_2_max_sb_herb = apply(pop_mats$seedbank_herb, MARGIN = 2, FUN = function(x) which(x == max(x))[1])

  df_out = as.data.frame(t(pop_mats$parameters))
  df_out$sb_over_yield = sb_over_yield
  df_out$time_2_oy_sb = as.numeric(time_2_oy_sb)
  df_out$time_in_oy_sb = as.numeric(time_in_oy_sb)
  df_out$mag_oy_sb = as.numeric(mag_oy_sb)
  df_out$min_diff_sb = as.numeric(diff_sb)
  df_out$fin_sb_herb = as.numeric(fin_sb_herb)
  df_out$fin_sb_noherb = as.numeric(fin_sb_noherb)
  df_out$max_sb_herb = as.numeric(max_sb_herb)
  df_out$max_sb_noherb = as.numeric(max_sb_noherb)
  df_out$time_2_max_sb_herb = as.numeric(time_2_max_sb_herb)
  df_out$time_2_max_sb_noherb = as.numeric(time_2_max_ag_noherb)
  df_out$ag_over_yield = ag_over_yield
  df_out$time_2_oy_ag = as.numeric(time_2_oy_ag)
  df_out$time_in_oy_ag = as.numeric(time_in_oy_ag)
  df_out$mag_oy_ag = as.numeric(mag_oy_ag)
  df_out$min_diff_ag = as.numeric(diff_ag)
  df_out$fin_ag_herb = as.numeric(fin_ag_herb)
  df_out$fin_ag_noherb = as.numeric(fin_ag_noherb)
  df_out$max_ag_herb = as.numeric(max_ag_herb)
  df_out$max_ag_noherb = as.numeric(max_ag_noherb)
  df_out$time_2_max_ag_herb = as.numeric(time_2_max_ag_herb)
  df_out$time_2_max_ag_noherb = as.numeric(time_2_max_ag_noherb)
  
  return(df_out)
}

## simplfied set of pop metrics used for filtering 
pop_metrics_simple <- function(pop_mats){

  df_out = as.data.frame(t(pop_mats$parameters))
  df_out$fin_sb_herb = as.numeric(pop_mats$seedbank_herb[dim(pop_mats$seedbank_herb)[1], ])
  df_out$fin_sb_noherb = as.numeric(pop_mats$seedbank_noherb[dim(pop_mats$seedbank_noherb)[1], ])
  df_out$fin_ag_noherb = as.numeric(pop_mats$ag_noherb[dim(pop_mats$ag_noherb)[1], ])
  df_out$fin_ag_pre_herb = as.numeric(pop_mats$ag_pre_herb[dim(pop_mats$ag_pre_herb)[1], ])
  df_out$fin_ag_post_herb = as.numeric(pop_mats$ag_post_herb[dim(pop_mats$ag_post_herb)[1], ])
  df_out$max_ag_post_herb = as.numeric(apply(pop_mats$ag_post_herb, MARGIN = 2, FUN = max))
  
  return(df_out)
}







