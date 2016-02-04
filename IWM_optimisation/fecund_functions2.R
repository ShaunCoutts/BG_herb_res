#you can ignore this function it is only there to setup the evaluation points which are a bit finicky, also will only be called onece at the start so has very little impact on the speed
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

# These next two functions need to be as fast as possible, they will be called many thousands of times and are also the slowest of the function used due to the implicit double loop
# at lines 40 (sum over first axis) and 64 (sum the result of that first sum, which is a 2D array, along a second axis to give a 1D array)
# These functions take a maternal and paternal distrbutions of individuals over a resistance score g (along with a bunch of parameters that control the proccess). 
# quant_gen_offspring_distribution() produces a 2D kernel (1 distribution over g of seeds produced for every value of g in the maternal distribution)
# This 2D distrbution is takn by fecundity() and reduced to a 1D kernel of seeds produced over evaluated values of g that includes all the varation from both the maternal and paternal 
# distrbutions of parents

## QUANT_GEN_OFFSPRING_DISTRIBUTION(N_f, eval_points, offspring_sd, seed_eval_points, dg) 
## produces a matrix where the rows are the distribution of offspring for each evaluated maternal breeding value based on paternal distributions on the breeeding vlaue (N_f) and 
##a vector of evaluation points on breeding value (each row in returned matrix is the distribution of offspring breeding values returned by each maternal breeding value evaluated), 
## along with parameter for the variance of breeding value for distribtuion of offspring breeding value (which we assume is normal).
## and a resolution to evluate the conditional offspring distribution at
#be careful this functions relies on N_f, being evaluated on eval_points before being passed to this function so that the indexes all match up
# Replaces original definition for easy benchmarking
quant_gen_offspring_distribution <- function(N_f, eval_points, offspring_sd, seed_eval_points, dg){
  eval_grid = cbind(rep.int(eval_points, length(eval_points)), rep(eval_points, each = length(eval_points)))#make every combination of evaluation points
  N_fathers = N_f / sum(N_f) #turns the distribution of fathers into a frequency distribution
  vect_seed_eval_points = rep(seed_eval_points, times = length(eval_grid[,1]))
  vect_breed_val_means = rep(eval_grid[,1] * 0.5 + eval_grid[,2] * 0.5, each = length(seed_eval_points))
  #MODIFIED: explicit calculation of vector of normal distributions
  vect_dnorm <- exp(-((vect_seed_eval_points - vect_breed_val_means)^2/(2*offspring_sd*offspring_sd)))/(sqrt(2*pi)*offspring_sd)
  cond_offspring_dist = matrix(vect_dnorm, ncol = length(seed_eval_points), byrow = TRUE) #creates matrix of normal distrubutions on each row 
  cond_offspring_dist = cond_offspring_dist * dg #scales the distribtuion so that it sums to one, the assumption being that all offspring produced by a parental combination have to have a breeding value
  offspring_3D_kernel = cond_offspring_dist * N_fathers #3D kernel held as a 2D matrix to make it easier and faster to work with 
  summing_grid = matrix(1:(length(eval_points) * length(eval_points)), ncol = length(eval_points), byrow = TRUE) #each row is a vector of indiceis that the conditional_offsprin_dist needs to be summed over to collapse the 3D kernel to 2D  
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
