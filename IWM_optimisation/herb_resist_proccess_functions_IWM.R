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

## QUANT_GEN_OFFSPRING_DISTRIBUTION(N_f, eval_points, additive_variance, offspring_dist_res) 
## produces a matrix where the rows are the distribution of offspring for each evaluated maternal breeding value based on paternal distributions on the breeeding vlaue (N_f) and 
##a vector of evaluation points on breeding value (each row in returned matrix is the distribution of offspring breeding values returned by each maternal breeding value evaluated), 
## along with parameter for the variance of breeding value for distribtuion of offspring breeding value (which we assume is normal).
## and a resolution to evluate the conditional offspring distribution at
#be careful this functions relies on N_f, being evaluated on eval_points before being passed to this function so that the indexes all match up
quant_gen_offspring_distribution <- function(N_f, eval_points, offspring_sd, seed_eval_points, dg){
  eval_grid = cbind(rep.int(eval_points, length(eval_points)), rep(eval_points, each = length(eval_points)))#make every combination of evaluation points
  N_fathers = N_f / sum(N_f) #turns the distribution of fathers into a frequency distribution
  vect_seed_eval_points = rep(seed_eval_points, times = length(eval_grid[,1]))
  vect_breed_val_means = rep(eval_grid[,1] * 0.5 + eval_grid[,2] * 0.5, each = length(seed_eval_points))
  cond_offspring_dist = matrix(dnorm(vect_seed_eval_points, vect_breed_val_means, offspring_sd), ncol = length(seed_eval_points), byrow = TRUE)
  cond_offspring_dist = cond_offspring_dist * dg #scales the distribtuion so that it sums to one, the assumption being that all offspring produced by a parental combination have to have a breeding value
  offspring_3D_kernel = cond_offspring_dist * N_fathers
  summing_grid = matrix(1:(length(eval_points) * length(eval_points)), ncol = length(eval_points), byrow = TRUE)
  t(apply(summing_grid, MARGIN = 1, FUN = function(x) colSums(offspring_3D_kernel[x, ])))
}

##FECUNDITY(N_m, eval_points, fec_max, fec0, fec_cost, N_f, additive_variance, seed_eval_points)
## produces a distribution of seeds on eval_points of g produced by population N_m (distribution of mothers on g evaluated at eval_points)
## N_m = maternal distrbution of indviduals over g
## N_f = paternal distrbution of indviduals over g, in most cases N_m == N_f
## eval_points = the values of g on which above ground individuals are evaluated
## seed_eval_points = the values of g on which seeds are evaluated
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## offspring_sd (passed to quant_gen_offspring_distribution()) = sd of conditional offspring distribution
## crop_effect_fec = proprional reduction in fecundity under different crops 
## density_effect_fec = reduction in density due to increased planting density 
## dense_depend_fec = 1/number of plants at which indivduals start to interfer with each other.
## dg = integration constant, equal to  eval_points[2] - eval_points[1]
#be careful this functions relies on N_m being evaluated on eval_points before being passed to this function so that the indexes all match up
fecundity <- function(N_m, eval_points, fec_max, fec0, fec_cost, N_f, offspring_sd, seed_eval_points, dense_depend_fec, crop_effect_fec, density_effect_fec, dg){
  num_survivors = sum(N_m) * dg
  if(num_survivors > 0){
    resist_effect_fec = exp(-(fec0 - eval_points * fec_cost))
    seeds_each_g = N_m * ((density_effect_fec * crop_effect_fec * fec_max) / (1 + resist_effect_fec + dense_depend_fec * num_survivors + 
      dense_depend_fec * num_survivors * resist_effect_fec))
    return(colSums(seeds_each_g * quant_gen_offspring_distribution(N_f, eval_points, offspring_sd, seed_eval_points, dg)) * dg )
  }else{
    return(0)
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
## mech_control = proportuion surviving mechanical control  
## crop_effect_sur = proportion surviving in a given grop relative to normal survival 
## pro_exposed = proprtion of emerged individuals that are exposed to herbicide
## seed_survival = probability that a seed in the seed bank survives one year
## seed_movement = movement of seed out of the seedbank in response to plowing
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## offspring_sd (passed to quant_gen_offspring_distribution()) = sd of conditional offspring distribution
## dense_depend_fec = density dependent effect on fecundity in units of 1/number individuals at which they start to affect each others fecundity
## crop_effect_fec = proprtion of normal fecundity in the crop choice choosen  
## density_effect_fec = proprtion of normal fecundity in the crop planting density choosen
## herb_rate = 0 or 1 factor that if herbicide was applied or not
## sur0 = susrvial rate when g is 0 and there is no herbicide
## sur_cost_resist = cost of higher resistenace score in terms of reduced survival
## herb_effect = effect of herbicide on survival
## survive_resist = protective effect of a one unit increase in resistance score g
## max_sur = maximum survival possible
## density_cutoff = density of emmerged plants at which are ignored in calculating the evaluation window, should be very small number
## dg = difference between evaluation points, should = eval_points[2] - eval_points[1]

single_iteration_1level <- function(seedbank_current, germination, mech_control, crop_effect_sur, seed_survival, seed_movement, eval_object, pro_exposed, herb_rate, sur0, 
    sur_cost_resist, herb_effect, survive_resist, max_sur, fec_max, fec0, fec_cost, offspring_sd, dense_depend_fec, density_cutoff, crop_effect_fec, 
    density_effect_fec, dg){
  
  new_seedbank = (seedbank_current - seedbank_current *  seed_movement) * seed_survival 
  new_plants = new_seedbank * germination * mech_control * crop_effect_sur   
  seedbank_post_germ = new_seedbank * (1 - germination)
  eval_object = eval_points_update(eval_points_object = eval_object, new_plants, density_cutoff = density_cutoff) #update evaluation window
  survivors_herb = pro_exposed * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground - herb_rate * (herb_effect - 
    pmin(herb_effect, survive_resist * eval_object$above_ground))))))
  survivors_noherb = (1 - pro_exposed) * new_plants[eval_object$above_ground_index] *  (max_sur / (1 + exp(-(sur0 - sur_cost_resist * eval_object$above_ground))))
  survivors_joint = survivors_herb + survivors_noherb
  seedbank_next = seedbank_post_germ + fecundity(N_m = survivors_joint, eval_points = eval_object$above_ground, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, 
    N_f = survivors_joint, offspring_sd = offspring_sd, seed_eval_points = eval_object$seed, dense_depend_fec = dense_depend_fec, crop_effect_fec = crop_effect_fec, 
    density_effect_fec = density_effect_fec, dg = dg) 
  
  return(seedbank_next)
}

## MULTI_ITERATION()
## produces a num_iter by length(eval_object$seed) matrix where each row is a distrbution over g of indivduals in the seedbank for 1 iteration   
## seedbank_current = distribution of seeds in the seed bank over g (either an intial seedbank or returned from another call to SINGLE_INTERATION()) 
## germination = germination probability
## mech_control = proportuion surviving mechanical control  
## crop_effect_sur = proportion surviving in a given grop relative to normal survival 
## pro_exposed = proprtion of emerged individuals that are exposed to herbicide
## seed_survival = probability that a seed in the seed bank survives one year
## seed_movement = movement of seed out of the seedbank in response to plowing
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## offspring_sd (passed to quant_gen_offspring_distribution()) = sd of conditional offspring distribution
## dense_depend_fec = density dependent effect on fecundity in units of 1/number individuals at which they start to affect each others fecundity
## crop_effect_fec = proprtion of normal fecundity in the crop choice choosen  
## density_effect_fec = proprtion of normal fecundity in the crop planting density choosen
## herb_rate = 0 or 1 factor that if herbicide was applied or not
## sur0 = susrvial rate when g is 0 and there is no herbicide
## sur_cost_resist = cost of higher resistenace score in terms of reduced survival
## herb_effect = effect of herbicide on survival
## survive_resist = protective effect of a one unit increase in resistance score g
## max_sur = maximum survival possible
## density_cutoff = density of emmerged plants at which are ignored in calculating the evaluation window, should be very small number
## dg = difference between evaluation points, should = eval_points[2] - eval_points[1]
## density_cutoff = population density (on distrbution over g) above which seed evaluation points are retained in the above ground evaluation points
## num_iter = number of interations to run the model for
## sub_action = an object the encodes all the actions into parameter values
## action_seq = matrix 5 x num_iter matrix that gives the sub-action taken in every iteration, used to reference action_space
multi_iteration <- function(seedbank_initial, germination, seed_survival, eval_object, pro_exposed, sur0, sur_cost_resist, herb_effect, survive_resist, max_sur, fec_max, 
  fec0, fec_cost, offspring_sd, dense_depend_fec, density_cutoff, dg, num_iter, sub_action, action_seq){
    
  results = matrix(NA, nrow = num_iter, ncol = length(eval_object$seed))
  i = 1
  results[i, ] = single_iteration_1level(seedbank_current = seedbank_initial, germination = germination, mech_control = sub_action$mech_sur[action_seq[i, 'mech']], 
    crop_effect_sur = sub_action$crop_sur[action_seq[i, 'crop']], seed_survival = seed_survival, seed_movement = sub_action$plow[action_seq[i, 'plow']], eval_object = eval_object, 
    pro_exposed = pro_exposed, herb_rate = sub_action$herb[action_seq[i, 'herb']], sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, 
    survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, 
    density_cutoff = density_cutoff, crop_effect_fec = sub_action$crop_fec[action_seq[i, 'crop']], density_effect_fec = sub_action$dens_fec[action_seq[i, 'dens']], dg = dg)

  for(i in 2:num_iter){
    results[i, ] = single_iteration_1level(seedbank_current = results[i - 1, ], germination = germination, mech_control = sub_action$mech_sur[action_seq[i, 'mech']], 
      crop_effect_sur = sub_action$crop_sur[action_seq[i, 'crop']], seed_survival = seed_survival, seed_movement = sub_action$plow[action_seq[i, 'plow']], eval_object = eval_object, 
      pro_exposed = pro_exposed, herb_rate = sub_action$herb[action_seq[i, 'herb']], sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, 
      survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, 
      density_cutoff = density_cutoff, crop_effect_fec = sub_action$crop_fec[action_seq[i, 'crop']], density_effect_fec = sub_action$dens_fec[action_seq[i, 'dens']], dg = dg)
  }
  
  results
}

## redo the answer key for the reference tests as changed all the bits around and removed some functions
# nonspace_test_answer_key[[1]]$question = paste0('eval_all = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = 0.5, seed_expantion = 4)\n',
#   'N_f = dnorm(eval_all$above_ground, 0, 2)\neval_points = eval_all$above_ground\noffspring_sd = 0.7\nseed_eval_points = eval_all$seed\ndg = 0.5\n')
# eval(parse(text = nonspace_test_answer_key[[1]]$question))
# nonspace_test_answer_key[[1]]$answer = quant_gen_offspring_distribution(N_f = N_f, eval_points = eval_points, offspring_sd = offspring_sd, seed_eval_points = seed_eval_points, dg = dg)
# 
# nonspace_test_answer_key[[2]]$question = paste0('lower_eval_point = -10\nupper_eval_point = 10\nresolution = 0.5\nseed_expantion = 4\n')
# eval(parse(text = nonspace_test_answer_key[[2]]$question)) 
# nonspace_test_answer_key[[2]]$answer = eval_points_builder(lower_eval_point = lower_eval_point, upper_eval_point = upper_eval_point, resolution = resolution, seed_expantion = seed_expantion)
# 
# nonspace_test_answer_key[[3]]$question = paste0('dg = 0.5\neval_all = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = dg, seed_expantion = 4)\n',
#   'N_m = N_f = dnorm(eval_all$above_ground, 0, 2)\neval_points = eval_all$above_ground\nfec_max = 200\nfec0 = 0\nfec_cost = 0.1\noffspring_sd = 0.7\nseed_eval_points = eval_all$seed\n',
#   'dense_depend_fec = 0.002\ncrop_effect_fec = 1\ndensity_effect_fec = 1\n')
# eval(parse(text = nonspace_test_answer_key[[3]]$question)) 
# nonspace_test_answer_key[[3]]$answer = fecundity(N_m = N_m, eval_points = eval_points, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, N_f = N_f, offspring_sd = offspring_sd, seed_eval_points = seed_eval_points, 
#   dense_depend_fec = dense_depend_fec, crop_effect_fec = crop_effect_fec, density_effect_fec = density_effect_fec, dg = dg)
# 
# nonspace_test_answer_key[[4]]$question = paste0('dg = 0.5\neval_object = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = dg, seed_expantion = 3)\n',
#   'seedbank_current = 200 * dnorm(eval_object$seed, 2, 3)\ngermination = 0.8\nmech_control = 1\ncrop_effect_sur = 1\npro_exposed = 0.7\nseed_survival = 0.5\nseed_movement = 0.3\n',
#   'fec_max = 200\nfec0 = 0\nfec_cost = 0.1\noffspring_sd = 0.7\ndense_depend_fec = 0.002\ncrop_effect_fec = 1\ndensity_effect_fec = 1\nherb_rate = 1\nsur0 = 10\nsur_cost_resist = 0\n',
#   'herb_effect = 10\nsurvive_resist = 0.5\nmax_sur = 0.95\ndensity_cutoff = 0.00001\n')
# eval(parse(text = nonspace_test_answer_key[[4]]$question)) 
# nonspace_test_answer_key[[4]]$answer = single_iteration_1level(seedbank_current = seedbank_current, germination = germination, mech_control = mech_control, crop_effect_sur = crop_effect_sur, seed_survival = seed_survival, 
#   seed_movement = seed_movement, eval_object = eval_object, pro_exposed = pro_exposed, herb_rate = herb_rate, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, 
#   survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, 
#   density_cutoff = density_cutoff, crop_effect_fec = crop_effect_fec, density_effect_fec = dense_depend_fec, dg = dg)
#   
# nonspace_test_answer_key[[5]]$question = 'eval_points_object = nonspace_test_answer_key[[2]]$answer\nabove_ground_dist = dnorm(eval_points_object$seed, 1, 3)\ndensity_cutoff = 0.00001'  
# eval(parse(text = nonspace_test_answer_key[[5]]$question))
# nonspace_test_answer_key[[5]]$answer = eval_points_update(eval_points_object = eval_points_object, above_ground_dist = above_ground_dist, density_cutoff = density_cutoff)
#   
# nonspace_test_answer_key[[6]]$question = paste0('dg = 0.5\neval_object = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = dg, seed_expantion = 3)\n',
#   'seedbank_initial = 200 * dnorm(eval_object$seed, 2, 3)\ngermination = 0.8\npro_exposed = 0.7\nseed_survival = 0.5\nfec_max = 200\nfec0 = 0\nfec_cost = 0.1\noffspring_sd = 0.7\n',
#   'dense_depend_fec = 0.002\nsur0 = 10\nsur_cost_resist = 0\nherb_effect = 10\nsurvive_resist = 0.5\nmax_sur = 0.95\ndensity_cutoff = 0.00001\ncrop_sur_alt = 0.8\ncrop_sur_fal = 0\n',
#   'crop_fec_alt = 0.95\ndens_fec_effect = 0.9\nmech_sur_effect = 0.1\nplow_effect = 0.3\n',
#   'sub_action = list(herb = c(0, 1), crop_sur = c(1, crop_sur_alt, crop_sur_fal), mech_sur = c(1, mech_sur_effect), dens_fec = c(1, dens_fec_effect), crop_fec = c(1, crop_fec_alt, 1), plow = c(0, plow_effect))\n',
#   'num_iter = 10\naction_seq = cbind(herb = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2), crop = c(1, 1, 1, 2, 2, 3, 1, 1, 1, 1), mech = c(1, 2, 1, 1, 1, 1, 1, 1, 1, 2), dens = c(1, 1, 2, 1, 1, 1, 2, 2, 2, 1), plow = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2))\n')
# eval(parse(text = nonspace_test_answer_key[[6]]$question))
# nonspace_test_answer_key[[6]]$answer = multi_iteration(seedbank_initial = seedbank_initial, germination = germination, seed_survival = seed_survival, eval_object = eval_object, pro_exposed = pro_exposed, 
#   sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, fec_max = fec_max, fec0 = fec0, 
#   fec_cost = fec_cost, offspring_sd = offspring_sd, dense_depend_fec = dense_depend_fec, density_cutoff = density_cutoff, dg = dg, num_iter = num_iter, sub_action = sub_action, 
#   action_seq = action_seq)
#   
# save(nonspace_test_answer_key, file = 'nonspatial_model_test_answer_key.Rdata')
#   
  
  
  

  ##area to test things with throw away code#######################################################################################################################
#library(microbenchmark)

#speed_test = microbenchmark(quant_gen_offspring_distribution(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed),
#  quant_gen_offspring_distribution_vect(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed),
#  times = 100)
#speed_test #turns out the fully vectorised version is much a bit slower which was unexpected but possibly due to large number of multiplications and additions required

#out1 = quant_gen_offspring_distribution(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed)
#out2 = quant_gen_offspring_distribution_vect(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed)
#all.equal(out1, out2)



























