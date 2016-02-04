library(microbenchmark)
setwd("~/Projects/code-review/20160128_shaun/")
source("fecund_functions.R")

#running the function and testing the speed
low_res = 0.5
high_res = 0.2
eval_points_low_res = eval_points_builder(-10, 10, resolution = low_res, 3) 
eval_points_high_res = eval_points_builder(-10, 10, resolution = high_res, 3)
dist_parents_low = 200 * dnorm(eval_points_low_res$above_ground, 0, 2)
dist_parents_high = 200 * dnorm(eval_points_high_res$above_ground, 0, 2)


speed_test = microbenchmark(fecundity(N_m = dist_parents_low, eval_points = eval_points_low_res$above_ground, fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_low, 
    offspring_sd = 0.7, seed_eval_points = eval_points_low_res$seed, dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = low_res),
  fecundity(N_m = dist_parents_high, eval_points = eval_points_high_res$above_ground, fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_high, 
    offspring_sd = 0.7, seed_eval_points = eval_points_high_res$seed, dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = high_res),
  times = 100)
speed_test

source("fecund_functions2.R")

speed_test = microbenchmark(fecundity(N_m = dist_parents_low, eval_points = eval_points_low_res$above_ground, fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_low, 
    offspring_sd = 0.7, seed_eval_points = eval_points_low_res$seed, dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = low_res),
  fecundity(N_m = dist_parents_high, eval_points = eval_points_high_res$above_ground, fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_high, 
    offspring_sd = 0.7, seed_eval_points = eval_points_high_res$seed, dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = high_res),
  times = 100)
speed_test

source("fecund_functions3.R")

speed_test = microbenchmark(fecundity(N_m = dist_parents_low, eval_points = eval_points_low_res$above_ground, fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_low, 
    offspring_sd = 0.7, seed_eval_points = eval_points_low_res$seed, dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = low_res),
  fecundity(N_m = dist_parents_high, eval_points = eval_points_high_res$above_ground, fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_high, 
    offspring_sd = 0.7, seed_eval_points = eval_points_high_res$seed, dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = high_res),
  times = 100)
speed_test

