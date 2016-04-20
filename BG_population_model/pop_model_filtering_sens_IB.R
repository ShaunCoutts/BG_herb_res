## ICEBERG VERSION

# Run the population model and to see how the population behaves and if the populaiton 
# can over yeild. We also use some field data to filter out parameter combinations 
# that produce unrealistic outcomes. We then use a sensitivtiy anaylsis to determine 
# which parameters are important in driving population growth, time to K and amount of 
# over-yeilding.


#file locations
setwd('/Users/shauncoutts/BG_pop_model')
#libraries and source codes
library(lhs)
library(parallel)

#set of functions for herbicide resitance model covering various processes like genotype production and seed bank dynamics

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

## GET_MEAN_SD()
## takes a distribution and gets the mean, sd and total number of seeds in the seedbank
get_mean_sd <- function(dist, eval_points, dg){
  total_sum = sum(dist) * dg
  approx_mean = sum((dist / total_sum) * eval_points * dg)
  approx_sd = sqrt(sum(((eval_points - approx_mean)^2) * (dist / total_sum) * dg))
  return(list(approx_mean = approx_mean, approx_sd = approx_sd, total_pop = total_sum))
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


# population model parameter est and ranges
offspring_sd = 1
sur0 = 10

param_ranges = rbind(seed_sur = c(0.22, 0.79), germ_prob = c(0.45, 0.6), fec_max = c(30, 300), 
  fec_dd = c(0.0000001, 0.0006), fec0 = c(5, 10), pro_exposed = c(0.6, 1), fec_cost_mult = c(0, 2),
  herb_effect_mult = c(2, 4), sur_protect_mult = c(0.001, 2))

# model run setup 
intial_pop_size = 100
dg = 0.5
eval_object = eval_points_builder(lower_eval_point = -20, upper_eval_point = 20, resolution = dg, 
  seed_expantion = 3)
seedbank_initial = dnorm(eval_object$seed, 0, 1.391607) * intial_pop_size
num_iter = 100
num_par_comb = 20000

#generate a set of parameters evenly sampled over the parameter space
sampled_points = improvedLHS(n = num_par_comb, k = 9, dup = 3)
sampled_param = lapply(1:dim(sampled_points)[1], FUN = function(x){
  ((param_ranges[, 2] - param_ranges[, 1]) * sampled_points[x, ]) + param_ranges[, 1]
}) 

t_par <- system.time({
  parameter_space_test = mclapply(sampled_param, FUN = function(x){ 
    pop_run(param_vect = x, int_pop_size = intial_pop_size, int_sd_g = 1.4, sur0 = 10, 
      offspring_sd = 1, burnin_time = 100, num_iter = 50, eval_object = eval_object, dg = dg)
  }, mc.cores = 8, mc.allow.recursive = FALSE)
})

t_par

save(parameter_space_test, file = 'parameter_space_test_ZHUMAC.Rdata')




