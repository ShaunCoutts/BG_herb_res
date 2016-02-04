library(lineprof)
setwd("~/Projects/code-review/20160128_shaun/")
source("fecund_functions.R")

#running the function and testing the speed
low_res = 0.5
high_res = 0.2
eval_points_low_res = eval_points_builder(-10, 10, resolution = low_res, 3) 
eval_points_high_res = eval_points_builder(-10, 10, resolution = high_res, 3)
dist_parents_low = 200 * dnorm(eval_points_low_res$above_ground, 0, 2)
dist_parents_high = 200 * dnorm(eval_points_high_res$above_ground, 0, 2)

l1 <- lineprof(fecundity(N_m = dist_parents_high,
                         eval_points = eval_points_high_res$above_ground,
                         fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_high, 
                         offspring_sd = 0.7, seed_eval_points = eval_points_high_res$seed,
                         dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = high_res))
shine(l1)

## this suggests that most time is being spent on the line that sets up conditional offspring distribution
## cond_offspring_dist = matrix(dnorm(vect_seed_eval_points, vect_breed_val_means, offspring_sd), ncol = length(seed_eval_points), byrow = TRUE)


source("fecund_functions2.R")

l2 <- lineprof(fecundity(N_m = dist_parents_high,
                         eval_points = eval_points_high_res$above_ground,
                         fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_high, 
                         offspring_sd = 0.7, seed_eval_points = eval_points_high_res$seed,
                         dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = high_res))
shine(l2)

## This shows good improvement from explicitly calculating normal distribution

source("fecund_functions3.R")

l3 <- lineprof(fecundity(N_m = dist_parents_high,
                         eval_points = eval_points_high_res$above_ground,
                         fec_max = 100, fec0 = 0, fec_cost = 1, N_f = dist_parents_high, 
                         offspring_sd = 0.7, seed_eval_points = eval_points_high_res$seed,
                         dense_depend_fec = 0.002, crop_effect_fec = 1, density_effect_fec = 1, dg = high_res))
shine(l3)

## This isn't any better than the previous attempt
