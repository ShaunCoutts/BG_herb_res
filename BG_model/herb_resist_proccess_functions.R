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
## above_ground_dist = a distribution of above ground plants evaluated on seed eval points, probably prduced by a call to EMERGENCE() 
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
quant_gen_offspring_distribution <- function(N_f, eval_points, additive_variance, seed_eval_points){
  dN = eval_points[2] - eval_points[1]
  additive_sd = sqrt(additive_variance)
  eval_grid = cbind(rep.int(eval_points, length(eval_points)), rep(eval_points, each = length(eval_points)))#make every combination of evaluation points
  N_fathers = N_f / sum(N_f) #turns the distribution of fathers into a frequency distribution
  vect_seed_eval_points = rep(seed_eval_points, times = length(eval_grid[,1]))
  vect_breed_val_means = rep(eval_grid[,1] * 0.5 + eval_grid[,2] * 0.5, each = length(seed_eval_points))
  cond_offspring_dist = matrix(dnorm(vect_seed_eval_points, vect_breed_val_means, additive_sd), ncol = length(seed_eval_points), byrow = TRUE)
  cond_offspring_dist = cond_offspring_dist * dN #scales the distribtuion so that it sums to one, the assumption being that all offspring produced by a parental combination have to have a breeding value
  offspring_3D_kernel = cond_offspring_dist * N_fathers
  summing_grid = matrix(1:(length(eval_points) * length(eval_points)), ncol = length(eval_points), byrow = TRUE)
  t(apply(summing_grid, MARGIN = 1, FUN = function(x) colSums(offspring_3D_kernel[x, ])))
}

## SURVIVAL(N_0, eval_points, herb_rate, sur0, sur_cost_resist, herb_effect, survive_resist, max_sur, ceiling_pop)
## produces a distribution (i.e. vector on eval_points) of survivors after herbicide application from a distribution of indivduals that emerge from the seed bank (N_0). 
## eval_points = vector of g values to evaluate the populaiton over
## herb_rate = 0 or 1 factor that if herbicide was applied or not
## sur0 = susrvial rate when g is 0 and there is no herbicide
## sur_cost_resist = cost of higher resistenace score in terms of reduced survival
## herb_effect = effect of herbicide on survival
## survive_resist = protective effect of a one unit increase in resistance score g
## max_sur = maximum survival possible
#be careful this functions relies on N_0 being evaluated on eval_points before being passed to this function so that the indexes all match up
survival <- function(N_0, eval_points, herb_rate, sur0, sur_cost_resist, herb_effect, survive_resist, max_sur, ceiling_pop){
  plant_happiness_sur = sur0 - sur_cost_resist * eval_points - herb_rate * (herb_effect - pmin(herb_effect, survive_resist * eval_points)) #alt: herb_effect * exp(-survive_resist * abs(g) + g)
  density_independent_survival = max_sur / (1 + exp(-plant_happiness_sur))
  density_independent_establishment = N_0 * density_independent_survival
  if((sum(density_independent_establishment)) > ceiling_pop){
    N_1 = density_independent_establishment * (ceiling_pop / sum(density_independent_establishment))
  }else{
    N_1 = density_independent_establishment
  }
  return(N_1)
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
## additive_variance (passed to quant_gen_offspring_distribution()) = variance of conditional offspring distribution
#be careful this functions relies on N_m being evaluated on eval_points before being passed to this function so that the indexes all match up
fecundity <- function(N_m, eval_points, fec_max, fec0, fec_cost, N_f, additive_variance, seed_eval_points){
  dg = eval_points[2] - eval_points[1]
  plant_happiness_fec = fec0 - fec_cost * eval_points
  seeds_each_g = N_m * (fec_max / (1 + exp(-plant_happiness_fec)))
  colSums(seeds_each_g * quant_gen_offspring_distribution(N_f, eval_points, additive_variance, seed_eval_points)) * dg 
}
  
## SEEDBANK()
## produces a distribution of seeds in the seed bank over the eval_points on g. 
## seedbank0 = distrbution of seeds in the seedbank in the last timestep
## seed_survival = probability that a seed in the seed bank survives one year
## germination = the probability that a seed in the seedbank survies one timestep
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## ALL PASSED TO FECUNDITY()
## N_m = maternal distrbution of indviduals over g
## N_f = paternal distrbution of indviduals over g, in most cases N_m == N_f
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## additive_variance (passed to quant_gen_offspring_distribution()) = variance of conditional offspring distribution
seedbank <- function(seedbank0, seed_survival, germination, eval_object, N_m, fec_max, fec0, fec_cost, N_f, additive_variance){
  seedbank0 * seed_survival * (1 - germination) + fecundity(N_m = N_m, eval_points = eval_object$above_ground, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, N_f = N_f, additive_variance = additive_variance, seed_eval_points = eval_object$seed)
}

## EMERGENCE()
## produces a distribution of emerged indviduals  
## seedbank_current = distribution of seeds in the seed bank over g (returned from SEEDBANK()) 
## germination = germination probability
emergence <- function(seedbank_current, germination){
  seedbank_current * germination
}

## SINGLE_INTERATION()
## produces a distrbution over g of indivduals in the seed bank including survival, reproduction and emergence 
## seedbank_current = distribution of seeds in the seed bank over g (returned from SEEDBANK() or another call to SINGLE_INTERATION()) 
## germination = germination probability
## seed_survival = probability that a seed in the seed bank survives one year
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## additive_variance (passed to quant_gen_offspring_distribution()) = variance of conditional offspring distribution
## herb_rate = 0 or 1 factor that if herbicide was applied or not
## sur0 = susrvial rate when g is 0 and there is no herbicide
## sur_cost_resist = cost of higher resistenace score in terms of reduced survival
## herb_effect = effect of herbicide on survival
## survive_resist = protective effect of a one unit increase in resistance score g
## max_sur = maximum survival possible
## density_cutoff = population density (on distrbution over g) above which seed evaluation points are retained in the above ground evaluation points

single_iteration <- function(seedbank_current, germination, eval_object, herb_rate, sur0, sur_cost_resist, herb_effect, survive_resist, max_sur, ceiling_pop, seed_survival, fec_max, fec0, fec_cost, additive_variance, density_cutoff){
  new_plants = emergence(seedbank_current = seedbank_current, germination = germination) 
  eval_object = eval_points_update(eval_points_object = eval_object, new_plants, density_cutoff = density_cutoff) #update evaluation window
  survivors = survival(N_0 = new_plants[eval_object$above_ground_index], eval_points = eval_object$above_ground, herb_rate = herb_rate, sur0 = sur0, sur_cost_resist = sur_cost_resist, 
	  herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop) 
  new_seedbank = seedbank(seedbank0 = seedbank_current, seed_survival = seed_survival, germination = germination, eval_object = eval_object, 
	  N_m =  survivors, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, N_f =  survivors, additive_variance = additive_variance) 
  new_seedbank
}

## MULTI_ITERATION()
## produces a num_iter by length(eval_object$seed) matrix where each row is a distrbution over g of indivduals in the seedbank for 1 iteration   
## num_iter = number of iterations to run the simulation form
## initial_seedbank = distribution of seeds in the seed bank over g
## herb_schedual = vector of 0,1's of length num_iter that defines when herbicide is applied 
## germination = germination probability
## seed_survival = probability that a seed in the seed bank survives one year
## eval_object = object from EVAL_POINTS_BUILDER() that defines the above ground and below ground evaluation points
## fec_max = the maximum number of seeds per mothers
## fec0 = cost of resistance when g = 0, in logits
## fec_cost = reduction in fecundity each additional unit of g causes, in logits
## additive_variance (passed to quant_gen_offspring_distribution()) = variance of conditional offspring distribution
## sur0 = susrvial rate when g is 0 and there is no herbicide
## sur_cost_resist = cost of higher resistenace score in terms of reduced survival
## herb_effect = effect of herbicide on survival
## survive_resist = protective effect of a one unit increase in resistance score g
## max_sur = maximum survival possible
## density_cutoff = population density (on distrbution over g) above which seed evaluation points are retained in the above ground evaluation points
multi_iteration <- function(num_iter, initial_seedbank, herb_schedual, germination, eval_object, sur0, sur_cost_resist, herb_effect, survive_resist, max_sur, ceiling_pop, seed_survival, fec_max, fec0, fec_cost, additive_variance, density_cutoff){
  results = matrix(NA, nrow = num_iter, ncol = length(eval_object$seed))
  results[1, ] = single_iteration(seedbank_current = initial_seedbank, germination = germination, eval_object = eval_object, herb_rate = herb_schedual[1], sur0 = sur0, sur_cost_resist = sur_cost_resist,
    herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop, seed_survival = seed_survival, fec_max = fec_max, fec0 = fec0, 
    fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff) 
  for(i in 2:num_iter){
    results[i, ] = single_iteration(seedbank_current = results[i - 1, ], germination = germination, eval_object = eval_object, herb_rate = herb_schedual[i], sur0 = sur0, sur_cost_resist = sur_cost_resist,
      herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop, seed_survival = seed_survival, fec_max = fec_max, fec0 = fec0, 
      fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff)
  }
  
  results
}

seedbank_animator <- function(results_matrix, eval_object, herb_schedual, pause = 1, ...){
  max_value = max(results_matrix)
  for(i in 1:dim(results_matrix)[1]){
    plot(eval_object$seed, eval_object$seed, type = 'n', ylim = c(0, max_value), bty = 'n', xlab = 'resistance score', 
      main = paste0(ifelse(herb_schedual[i] == 0, 'No herbicide applied', 'Herbicide being applied'), '\nturn ', i), ...)
    polygon(x = eval_object$seed, y = results_matrix[i, ], col = ifelse(herb_schedual[i] == 0, 'skyblue', 'red'))
    Sys.sleep(pause)
  }
}


#Check the test results  
test_functions_broken <- function(file_loc){
  setwd(file_loc)
  test_obj_name <- load('nonspatial_model_test_answer_key.Rdata') #load the test key
  eval(parse(text = nonspace_test_answer_key[[1]]$question))#set parameters for the test run
  test1 = quant_gen_offspring_distribution(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed) #get the output form the current version of the function
  eval(parse(text = nonspace_test_answer_key[[2]]$question))#set parameters for the test run
  test2 = survival(N_0 = N_0, eval_points = eval_points, herb_rate = herb_rate, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop) 
  eval(parse(text = nonspace_test_answer_key[[3]]$question))#set parameters for the test run
  test3 = eval_points_builder(lower_eval_point = lower_eval_point, upper_eval_point = upper_eval_point, resolution = resolution, seed_expantion = seed_expantion)
  eval(parse(text = nonspace_test_answer_key[[4]]$question))#set parameters for the test run
  test4 = fecundity(N_m = N_m, eval_points = eval_all$above_ground, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, N_f = N_f, additive_variance = additive_variance, seed_eval_points = eval_all$seed) 
  eval(parse(text = nonspace_test_answer_key[[5]]$question))#set parameters for the test run
  test5 = seedbank(seedbank0 = seedbank0 * 100, seed_survival = seed_survival, germination = germination, eval_object = eval_object, N_m = N_m, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, N_f = N_f, additive_variance = additive_variance) 
  eval(parse(text = nonspace_test_answer_key[[6]]$question))#set parameters for the test run
  test6 = emergence(seedbank_current = seedbank_current, germination = germination) 
  eval(parse(text = nonspace_test_answer_key[[7]]$question))#set parameters for the test run
  test7 = single_iteration(seedbank_current = seedbank_current, germination = germination, eval_object = eval_object, herb_rate = herb_rate, sur0 = sur0, sur_cost_resist = sur_cost_resist,
    herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop, seed_survival = seed_survival, fec_max = fec_max, fec0 = fec0, 
    fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff) 
  eval(parse(text = nonspace_test_answer_key[[8]]$question))#set parameters for the test run
  test8 = eval_points_update(eval_points_object = eval_points_object, above_ground_dist = above_ground_dist, density_cutoff = density_cutoff)
  eval(parse(text = nonspace_test_answer_key[[9]]$question))#set parameters for the test run
  test9 = multi_iteration(num_iter = num_iter, initial_seedbank = initial_seedbank, germination = germination, eval_object = eval_object, herb_schedual = herb_schedual, sur0 = sur0, sur_cost_resist = sur_cost_resist,
    herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop, seed_survival = seed_survival, fec_max = fec_max, fec0 = fec0, 
    fec_cost = fec_cost, additive_variance = additive_variance, density_cutoff = density_cutoff) 

  test_results <- ifelse(all.equal(test1, nonspace_test_answer_key[[1]]$answer), 'QUANT_GEN_OFFSPRING_DISTRIBUTION() still fine', 'Something you did broke the function QUANT_GEN_OFFSPRING_DISTRIBUTION()')
  test_results[2] <- ifelse(identical(test2, nonspace_test_answer_key[[2]]$answer), 'SURVIVAL() still fine', 'Something you did broke the function SURVIVAL()')
  test_results[3] <- ifelse(identical(test3, nonspace_test_answer_key[[3]]$answer), 'EVAL_POINTS_BUILDER() still fine', 'Something you did broke the function EVAL_POINTS_BUILDER()')
  test_results[4] <- ifelse(identical(test4, nonspace_test_answer_key[[4]]$answer), 'FECUNDITY() still fine', 'Something you did broke the function FECUNDITY()')
  test_results[5] <- ifelse(identical(test5, nonspace_test_answer_key[[5]]$answer), 'SEEDBANK() still fine', 'Something you did broke the function SEEDBANK()')
  test_results[6] <- ifelse(identical(test6, nonspace_test_answer_key[[6]]$answer), 'EMERGENCE() still fine', 'Something you did broke the function EMERGENCE()')
  test_results[7] <- ifelse(identical(test7, nonspace_test_answer_key[[7]]$answer), 'SINGLE_INTERATION() still fine', 'Something you did broke the function SINGLE_INTERATION()')
  test_results[8] <- ifelse(identical(test8, nonspace_test_answer_key[[8]]$answer), 'EVAL_POINTS_UPDATE() still fine', 'Something you did broke the function EVAL_POINTS_UPDATE()')
  test_results[9] <- ifelse(identical(test9, nonspace_test_answer_key[[9]]$answer), 'MULTI_ITERATION() still fine', 'Something you did broke the function MULTI_ITERATION()')

  print(test_results) 
}





#NOTE TO SELF add non-heritable variance in resitance in the fecundity function so individuals can be resistant through their life time, basically n(g) should be n(g, z) and resistance should then 
#be a function of g and z, so that survival is actually a distribtuion for each element of eval_points do simple version for now.

  
  
  
  
  
  


  
  
  
  

  
  
  ##area to test things with throw away code#######################################################################################################################
#library(microbenchmark)

#speed_test = microbenchmark(quant_gen_offspring_distribution(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed),
#  quant_gen_offspring_distribution_vect(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed),
#  times = 100)
#speed_test #turns out the fully vectorised version is much a bit slower which was unexpected but possibly due to large number of multiplications and additions required

#out1 = quant_gen_offspring_distribution(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed)
#out2 = quant_gen_offspring_distribution_vect(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed)
#all.equal(out1, out2)



























