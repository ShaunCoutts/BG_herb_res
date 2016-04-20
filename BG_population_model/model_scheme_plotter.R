library(colorspace)

#get functions from the population model to help with plotting 
setwd('/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code')
source('herb_resist_proccess_functions_IWM.R')

fec_func = 
# function to produce a shematic plot of the population to show how population distribution,
# survival and fecudity functions relate to each other 
pop_scheme <- function(pop_sd, mean_g, eval_points_object, fec0, fec_cost, sur0, herb_rate, herb_effect, survive_resist, 
  col_pop, col_sur, col_fec, dg){
  
  fec_func = fecundity_closure(eval_points_object, offspring_sd = 1, fec0 = fec0, fec_cost = fec_cost)
  
  x_points = eval_points_object$above_ground
  par(mar = c(4.5, 0.5, 1.5, 0.5), mfrow = c(2, 1))
  plot(x_points, seq(0, 1.15, length = length(x_points)), type = 'n', bty = 'n', yaxt = 'n', ylab = '', ylim = c(0, 1.15), xlab = 'resistance')
  lines(x = c(0, 0), y = c(-0.2, 1.05))
  text(x = 0, y = 1.1, labels = 'prob. or\nfreq.', pos = 2)
  
  pop = dnorm(x_points, mean_g, pop_sd)	
  scale_fact = 1 / max(pop)
  pop = pop * scale_fact
  polygon(x = x_points, y = pop, col = adjustcolor(col_pop, alpha.f = 0.7), border = col_pop)
  text(x = x_points[floor(length(x_points) * 0.70)], y = pop[floor(length(pop) * 0.70)], labels = 'population', col = col_pop, pos = 4)
 
  #add survival curve 
  sur_curve = 1 / (1 + exp(-(sur0 - herb_rate * (herb_effect - pmin(herb_effect, survive_resist * x_points)))))
  lines(x_points, sur_curve, lwd = 2, col = col_sur)
  text(x = x_points[floor(length(x_points) * 0.70)], y = sur_curve[floor(length(sur_curve) * 0.70)] - 0.03, labels = 'survival', col = col_sur, pos = 4)
  
  #plot survivours and seeds produced from them
  plot(x_points, seq(0, 1.15, length = length(x_points)), type = 'n', bty = 'n', yaxt = 'n', ylab = '', ylim = c(0, 1.15), xlab = 'resistance')
  lines(x = c(0, 0), y = c(-0.2, 1.05))
  text(x = 0, y = 1.1, labels = 'prob. or\nfreq.', pos = 2)
  #add in seed distribution
  new_seeds = fec_func(N_m = pop_sur, fec_max = 1, dense_depend_fec = 0, density_effect_fec = 1, crop_effect_fec = 1, dg = dg)
  new_seeds = new_seeds[eval_points_object$above_ground_index] * (1 / max(new_seeds[eval_points_object$above_ground_index]))
  polygon(x = x_points, y = new_seeds, col = adjustcolor(col_fec, alpha.f = 0.7), border = col_fec)
  text(x = x_points[floor(length(x_points) * 0.70)], y = new_seeds[floor(length(new_seeds) * 0.70)], labels = 'seeds', col = col_fec, pos = 4)
  #add survivors distribution
  pop_sur = pop * sur_curve 
  polygon(x = x_points, y = pop_sur, col = adjustcolor(col_pop, alpha.f = 0.7), border = col_pop)
  text(x = x_points[floor(length(x_points) * 0.70)], y = pop_sur[floor(length(pop_sur) * 0.70)] + 0.05, labels = 'population\nsurvivors', col = col_pop, pos = 4)
  # add in fec line 
  fec_curve = 1 / (1 + exp(-(fec0 - fec_cost * x_points)))
  lines(x_points, fec_curve, lwd = 2, col = col_fec)
  text(x = x_points[floor(length(x_points) * 0.30)], y = fec_curve[floor(length(fec_curve) * 0.30)] - 0.03, labels = 'fecundity', col = col_fec, pos = 4)
  
  return(NA)
}

dg = 0.05    
eval_points_object = eval_points_builder(lower_eval_point= -6, upper_eval_point = 6, resolution = dg, seed_expantion = 2)
pop_sd = 1
mean_g = 0
col_pop = rainbow_hcl(n = 1, start = 140)
col_sur = rainbow_hcl(n = 1, start = 0)
col_fec = rainbow_hcl(n = 1, start = 270)
sur0 = 3
herb_rate = 1
herb_effect = 5
survive_resist = 4
fec0 = 4 
fec_cost = 5


pop_scheme(pop_sd = 1.4, mean_g = 0, eval_points_object = eval_points_object, fec0 = 4, fec_cost = 5, sur0 = 3, herb_rate = 1, herb_effect = 5, survive_resist = 4, 
  col_pop = col_pop, col_sur = col_sur, col_fec = col_fec, dg = dg){