#makes the evaluation mesh
eval_points_builder <- function(lower_eval_point, upper_eval_point, resolution, seed_expantion){
  above_ground_eval = seq(lower_eval_point, upper_eval_point, resolution)
  seed_lower = seq(above_ground_eval[1] * seed_expantion, above_ground_eval[1], resolution)
  seed_upper = seq(above_ground_eval[length(above_ground_eval)], above_ground_eval[length(above_ground_eval)] * seed_expantion, resolution)
  seed_eval = c(seed_lower[1:(length(seed_lower) - 1)], above_ground_eval, seed_upper[2:length(seed_upper)])
  list(above_ground = above_ground_eval, seed = seed_eval, above_ground_index = which(seed_eval %in% above_ground_eval))
}
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
#makes the quant_gen function to call later
dg = 0.5
eval_points_obj = eval_points_builder(lower_eval_point = -10, upper_eval_point = 10, resolution = dg, seed_expantion = 3)
offspring_sd = 0.7

quant_gen_offspring_distribution = quant_gen_closure(eval_points_object = eval_points_obj, offspring_sd = offspring_sd)
#makes the distribution of new seeds over g with cost of resistance and density and other control options affecting seed number 
fecundity <- function(N_m, eval_points, fec_max, fec0, fec_cost, dense_depend_fec, crop_effect_fec, density_effect_fec, dg){
  num_survivors = sum(N_m) * dg
  if(num_survivors > 0){
    resist_effect_fec = exp(-(fec0 - eval_points * fec_cost))
    repo_happiness = (density_effect_fec * crop_effect_fec) / (1 + resist_effect_fec + dense_depend_fec * num_survivors + dense_depend_fec * num_survivors * resist_effect_fec)
    post_select_parents = N_m * repo_happiness
    post_select_parents_normz = post_select_parents / (sum(post_select_parents) * dg)
    pd_seeds = quant_gen_offspring_distribution(N_m = post_select_parents_normz) * dg * dg * dg #need to integrate 3 times across mix_kern cols, N_m and N_p so whole thing sums to 1
    return(sum(repo_happiness * post_select_parents) * fec_max * dg * pd_seeds) #some slight inaccuracy here, if the expected number of total seeds is 20,000 the actual total number is 19,998.4
 }else{
    return(0)
  }
}


#parameters
library(microbenchmark)
N_f = 200 * dnorm(eval_points_obj$above_ground, 0, 0.98)
fec_max = 100 
fec0 = 100 
fec_cost = 0 
dense_depend_fec = 0.0000002 
crop_effect_fec = 1 
density_effect_fec = 1 

microbenchmark(fecundity(N_m = N_f, eval_points = eval_points_obj$above_ground, fec_max = fec_max, fec0 = fec0, fec_cost = fec_cost, dense_depend_fec = dense_depend_fec, 
  crop_effect_fec = crop_effect_fec, density_effect_fec = density_effect_fec, dg = dg), times = 20)
 