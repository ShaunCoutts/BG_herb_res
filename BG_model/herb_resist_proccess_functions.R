#set of functions for herbicide resitance model covering various processes like genotype production and seed bank dynamics
working_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_model' 
setwd(working_loc)
test_obj_name <- load('nonspatial_model_test_answer_key.Rdata') #load the test key

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
eval(parse(text = nonspace_test_answer_key[[3]]$question))#set parameters for the test run
test3 = eval_points_builder(lower_eval_point = lower_eval_point, upper_eval_point = upper_eval_point, resolution = resolution, seed_expantion = seed_expantion)
#nonspace_test_answer_key[[3]] = list(question = 'lower_eval_point = -10\nupper_eval_point = 10\nresolution = 0.5\nseed_expantion = 4\n', answer = test3)
#setwd(working_loc)
#save(nonspace_test_answer_key, file = 'nonspatial_model_test_answer_key.Rdata')

## QUANT_GEN_OFFSPRING_DISTRIBUTION(N_f, eval_points, additive_variance, offspring_dist_res) 
## produces a matrix where the rows are the distribution of offspring for each evaluated maternal breeding value based on paternal distributions on the breeeding vlaue (N_f) and 
##a vector of evaluation points on breeding value (each row in returned matrix is the distribution of offspring breeding values returned by each maternal breeding value evaluated), 
## along with parameter for the variance of breeding value for distribtuion of offspring breeding value (which we assume is normal).
## and a resolution to evluate the conditional offspring distribution at
#be careful this functions relies on N_f, being evaluated on eval_points before being passed to this function so that the indexes all match up
quant_gen_offspring_distribution <- function(N_f, eval_points, additive_variance, seed_eval_points){
  additive_sd = sqrt(additive_variance)
  eval_grid = cbind(rep.int(eval_points, length(eval_points)), rep(eval_points, each = length(eval_points)))#make every combination of evaluation points
  index_comb = cbind(rep.int(1:length(eval_points), length(eval_points)), rep(1:length(eval_points), each = length(eval_points)))#make every combination of index on the 1st and 2nd dimention of eval_grid
  N_fathers = N_f / sum(N_f) #turns the distribution of fathers into a frequency distribution
  offspring_3D_kernel = sapply(seq_along(index_comb[,1]), FUN = function(x){
    cond_offspring_dist = dnorm(seed_eval_points, 0.5 * eval_grid[x, 1] + 0.5 * eval_grid[x, 2], additive_sd) #centers the conditional distribtuion of offspring on the breeding value being assesed so that at the extreams the full distribution is still assesed
    cond_offspring_dist = cond_offspring_dist / sum(cond_offspring_dist) #scales the distribtuion so that it sums to one, the assumption being that all offspring produced by a parental combination have to have a breeding value
    cond_offspring_dist * N_fathers[index_comb[x, 1]]  
  })
  offspring_kernel_dims = dim(offspring_3D_kernel)
  summing_grid = matrix(1:(length(eval_points) * length(eval_points)), ncol = length(eval_points), byrow = TRUE)
  sapply(1:offspring_kernel_dims[1], FUN = function(i) apply(summing_grid, MARGIN = 1, FUN = function(x) sum(offspring_3D_kernel[i, x])))
}
#test of quant_gen_offspring_distribution() with dummy input and a fixed known output so I can see if I break this at some point in the future
eval(parse(text = nonspace_test_answer_key[[1]]$question))#set parameters for the test run
test1 = quant_gen_offspring_distribution(N_f = N_f, eval_points = eval_all$above_ground, additive_variance = additive_variance, seed_eval_points = eval_all$seed) #get the output form the current version of the function
#evaluate this aginst the reference answer
#write the output to a named list stored on disk so that this functions output can be tested aginast it at future date if it is changed.
#nonspace_test_answer_key = list()
#nonspace_test_answer_key[[1]] = list(question = 'eval_all = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = 1.5, seed_expantion = 4)\nN_f = dnorm(eval_all$above_ground, 0, 2)\nadditive_variance = 0.5\n', answer = test1)
#setwd(working_loc)
#save(nonspace_test_answer_key, file = 'nonspatial_model_test_answer_key.Rdata')

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
  plant_happiness_sur = sur0 - sur_cost_resist * eval_points - herb_effect * herb_rate + survive_resist * eval_points * herb_rate #linear function that combines cost of resistance, effect of herbicide and the protective effect of g
  density_independent_survival = max_sur / (1 + exp(-plant_happiness_sur))
  density_independent_establishment = N_0 * density_independent_survival
  ifelse(sum(density_independent_establishment) > ceiling_pop, density_independent_establishment * (ceiling_pop / sum(density_independent_establishment)), density_independent_establishment)
}

#evaluate this aginst the reference answer
eval(parse(text = nonspace_test_answer_key[[2]]$question))#set parameters for the test run
test2 = survival(N_0 = N_0, eval_points = eval_points, herb_rate = herb_rate, sur0 = sur0, sur_cost_resist = sur_cost_resist, herb_effect = herb_effect, survive_resist = survive_resist, max_sur = max_sur, ceiling_pop = ceiling_pop) 
#put the question and answer in the answer key
#nonspace_test_answer_key[[2]] = list(question = 'eval_points = seq(-10, 10, 1.5)\nN_0 = dnorm(eval_points, 0, 2)\nherb_rate = 1\nsur0 = 5\nsur_cost_resist = 0.1\nherb_effect = 3\nsurvive_resist = 5\nmax_sur = 0.95\nceiling_pop = 2\n', answer = test2)
#save(nonspace_test_answer_key, file = 'nonspatial_model_test_answer_key.Rdata')

##FECUNDITY()
## produces a distribution of seeds on eval_points of g produced by population N_m (distribution of mothers on g evaluated at eval_points) 
#be careful this functions relies on N_m being evaluated on eval_points before being passed to this function so that the indexes all match up
fecundity <- function(N_m, eval_points, fec_max, fec0, fec_cost, N_f, additive_variance, seed_eval_points){
  plant_happiness_fec = fec0 - fec_cost * eval_points
  seeds_each_g = N_m * (fec_max / (1 + exp(-plant_happiness_fec)))
  colSums(seeds_each_g * quant_gen_offspring_distribution(N_f, eval_points, additive_variance, seed_eval_points)) #this needs to be checked, not 100% sure
  
}


## SEEDBANK()
## produces a distribution of seeds in the seed bank over the eval_points on g. 



  
  
#NOTE TO SELF add non-heritable variance in resitance in the fecundity function so individuals can be resistant through their life time, basically n(g) should be n(g, z) and resistance should then 
#be a function of g and z, so that survival is actually a distribtuion for each element of eval_points do simple version for now.
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#Check the test results  
test_results <- ifelse(identical(test1, nonspace_test_answer_key[[1]]$answer), 'QUANT_GEN_OFFSPRING_DISTRIBUTION() still fine', 'Something you did broke the function QUANT_GEN_OFFSPRING_DISTRIBUTION()')
test_results[2] <- ifelse(identical(test2, nonspace_test_answer_key[[2]]$answer), 'SURVIVAL() still fine', 'Something you did broke the function SURVIVAL()')
test_results[3] <- ifelse(identical(test3, nonspace_test_answer_key[[3]]$answer), 'EVAL_POINTS_BUILDER() still fine', 'Something you did broke the function EVAL_POINTS_BUILDER()')
print(test_results) 
  
  
  
  

  
  
  ##area to test things with throw away code#######################################################################################################################
mat = cbind(g_m = rep.int(1:10, length(1:10)), g_f = rep(1:10, each = length(1:10)))
aperm('dim<-' (mat, list(10, 10, 2)), c(2, 1, 3))
