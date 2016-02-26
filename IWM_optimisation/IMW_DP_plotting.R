#plotting scripts and functions for dynamic programing output
library(shape)
#simple plot of best action in action space on the sampled state space
plot_policy <- function(Q, output_loc, output_name = 'policy_output.pdf'){
  inital_loc = getwd()
  setwd(output_loc)
  col_palett = rainbow_hcl(n = 25, start = 0, end = 200) 
  #set up an evaluation grid to plot the value surface from the gam 
  upper_mean = max(Q[[1]]$sample_points_df$mean_g)
  lower_mean = min(Q[[1]]$sample_points_df$mean_g)
  upper_pop = max(Q[[1]]$sample_points_df$pop)
  lower_pop = min(Q[[1]]$sample_points_df$pop)
  eval_res = 100
  value_grid = data.frame(mean_g = rep.int(seq(lower_mean, upper_mean, length.out = eval_res), eval_res), pop = rep(seq(lower_pop, upper_pop, length.out = eval_res), each = eval_res))#make every combination of evaluation points
  
  pdf(file = output_name, width = 15, height = 10)
    for(i in seq_along(Q)){
      par(mfrow = c(1, 3))
      plot(x = Q[[i]]$sample_points_df$mean_g, y = Q[[i]]$sample_points_df$pop, type = 'n', xlab = 'mean_g', ylab = 'population', main = paste0('Policy t = T - ', (i - 1)))
      text(x = Q[[i]]$sample_points_df$mean_g , y = Q[[i]]$sample_points_df$pop , labels = paste0('a', Q[[i]]$sample_points_df$policy), 
	col = col_palett[Q[[i]]$sample_points_df$policy])
      #plot the value surface
      if(class(Q[[i]]$value_surface) == 'gbm'){
	cat('plotting surface', i, '\n')
	image(matrix(predict(Q[[i]]$value_surface, newdata = value_grid, n.trees = Q[[i]]$value_surface$n.trees), ncol = eval_res, byrow = FALSE), xaxt = 'n', xlab = 'mean_g', yaxt = 'n', ylab = 'population')
      }else{
	image(matrix(predict(Q[[i]]$value_surface, value_grid), ncol = eval_res, byrow = FALSE), xaxt = 'n', xlab = 'mean_g', yaxt = 'n', ylab = 'population')
      }
      axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(value_grid$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1))
      axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(value_grid$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 1))
      #plot fitted VS predicted for spline
      if(class(Q[[i]]$value_surface) == 'gbm'){
	plot(x = Q[[i]]$sample_points_df$values, y = Q[[i]]$value_surface$fit ,xlab = 'observed', ylab = 'predicted') 
	abline(0, 1)
      }else{
	plot(x = Q[[i]]$sample_points_df$values, y = predict(Q[[i]]$value_surface, Q[[i]]$sample_points_df), xlab = 'observed', ylab = 'predicted') 
	abline(0, 1)
      }
    }
  dev.off()
  setwd(inital_loc)
}

#a much prettier plot of the policy with none of the fitting images included
#named_col_pal must be a vector of colors named with the class names, so I can keep colours consistent between plots 
pretty_policy_plot <- function(Q, res, named_col_pal, output_loc, output_name = 'policy_pretty_output.pdf', ...){
  mean_g_range = range(Q$sample_points_df$mean_g)
  pop_range = range(Q$sample_points_df$pop)
  #rescale so both state variables on 0 - 1 scale 
  mean_g_scaled = (Q$sample_points_df$mean_g - mean_g_range[1]) / (mean_g_range[2] - mean_g_range[1]) 
  pop_scaled = (Q$sample_points_df$pop - pop_range[1]) / (pop_range[2] - pop_range[1]) 
  #size for squares
  rect_size = 1 / res
  #set up layout matrix
  layout_mat = rbind(c(0, 0.75, 0, 1), c(0.75, 1, 0, 1))
  #make every combination of mean_g and pop in mesh with res = res 
  test_points = data.frame(mean_g = rep.int(seq(0, 1, length.out = res), res), pop = rep(seq(0, 1, length.out = res), each = res))#make every combination of evaluation points
  pred_class = as.character(knn1(train = data.frame(mean_g = mean_g_scaled, pop = pop_scaled), test = test_points, 
    cl = paste0('a_', Q$sample_points_df$policy)))
    
  current_loc = getwd()
  setwd(output_loc)
  pdf(file = output_name, width = 15, height = 11)
    split.screen(layout_mat)
    screen(1)
      par(mar = c(5, 5, 1, 1))
      plot(x = 0:1, y = 0:1, type = 'n', xlim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size), ylim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size),
	bty = 'n', tck = 0.02, xaxt = 'n', yaxt = 'n', ...)
      axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1), tck = 0.015, 
	cex.axis = 1.5)
      axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 1), tck = 0.015, 
	cex.axis = 1.5)
      
      rect(xleft = test_points$mean_g - 0.505 * rect_size, xright = test_points$mean_g + 0.505 * rect_size, ybottom = test_points$pop - 0.505 * rect_size, 
	ytop = test_points$pop + 0.505 * rect_size, border = NA, col = adjustcolor(named_col_pal[pred_class], alpha.f = 0.4))
      
      text(x = mean_g_scaled, y = pop_scaled, labels = paste0('a_', Q$sample_points_df$policy), col = named_col_pal[paste0('a_', Q$sample_points_df$policy)])
    screen(2)
      plot(x = 0:1, y = 0:1, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
      legend(x = 'center', legend = paste0('action ', unique(Q$sample_points_df$policy)), fill = named_col_pal[paste0('a_', unique(Q$sample_points_df$policy))],
	cex = 2, bty = 'n', border = NA)
    close.screen(all.screens = TRUE)
  dev.off()   
  setwd(current_loc)
}

#function to tease out the different sub actions in each policy with 1 plot for each sub-action
policy_breakdown_plot <- function(Q, action_space, res, breakdown_col_pal, output_loc, output_name = 'policy_breakdown_output.pdf', ...){
  mean_g_range = range(Q$sample_points_df$mean_g)
  pop_range = range(Q$sample_points_df$pop)
  #rescale so both state variables on 0 - 1 scale 
  mean_g_scaled = (Q$sample_points_df$mean_g - mean_g_range[1]) / (mean_g_range[2] - mean_g_range[1]) 
  pop_scaled = (Q$sample_points_df$pop - pop_range[1]) / (pop_range[2] - pop_range[1]) 
  #size for squares
  rect_size = 1 / res
  #make every combination of mean_g and pop in mesh with res = res 
  test_points = data.frame(mean_g = rep.int(seq(0, 1, length.out = res), res), pop = rep(seq(0, 1, length.out = res), each = res))#make every combination of evaluation points
  #make the prediction on the test points from the sample points
  pred_action = as.numeric(as.character(knn1(train = data.frame(mean_g = mean_g_scaled, pop = pop_scaled), test = test_points, 
    cl = Q$sample_points_df$policy)))
    
  pred_crop = action_space[pred_action, 'crop']
  pred_herb = action_space[pred_action, 'herb']
  pred_plow = action_space[pred_action, 'plow']
  pred_mech = action_space[pred_action, 'mech']
  pred_dens = action_space[pred_action, 'dens']
  
  current_loc = getwd()
  setwd(output_loc)
  pdf(file = output_name, width = 10.1, height = 15)
    par(mfcol = c(3, 2), mar = c(5, 5, 3, 1))
    #crop action
    plot(x = 0:1, y = 0:1, type = 'n', xlim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size), ylim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size),
      bty = 'n', tck = 0.02, xaxt = 'n', yaxt = 'n', main = 'Crop choice', cex.main = 2, cex.lab = 2, ...)
    axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1), tck = 0.015, 
      cex.axis = 1.5)
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 0), tck = 0.015, 
      cex.axis = 1.5)
    rect(xleft = test_points$mean_g - 0.505 * rect_size, xright = test_points$mean_g + 0.505 * rect_size, ybottom = test_points$pop - 0.505 * rect_size, 
      ytop = test_points$pop + 0.505 * rect_size, border = NA, col = adjustcolor(breakdown_col_pal$crop[pred_crop], alpha.f = 0.7))
    points(x = mean_g_scaled , y = pop_scaled, pch = 19, cex = 1.5, col = breakdown_col_pal$crop[action_space[Q$sample_points_df$policy, 'crop']]) 
    #herbicide 
    plot(x = 0:1, y = 0:1, type = 'n', xlim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size), ylim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size),
      bty = 'n', tck = 0.02, xaxt = 'n', yaxt = 'n', main = 'Herbicide', cex.main = 2, cex.lab = 2, ...)
    axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1), tck = 0.015, 
      cex.axis = 1.5)
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 0), tck = 0.015, 
      cex.axis = 1.5)
    rect(xleft = test_points$mean_g - 0.505 * rect_size, xright = test_points$mean_g + 0.505 * rect_size, ybottom = test_points$pop - 0.505 * rect_size, 
      ytop = test_points$pop + 0.505 * rect_size, border = NA, col = adjustcolor(breakdown_col_pal$herb[pred_herb], alpha.f = 0.7))
    points(x = mean_g_scaled , y = pop_scaled, pch = 19, cex = 1.5, col = breakdown_col_pal$herb[action_space[Q$sample_points_df$policy, 'herb']]) 
    #plow
    plot(x = 0:1, y = 0:1, type = 'n', xlim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size), ylim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size),
      bty = 'n', tck = 0.02, xaxt = 'n', yaxt = 'n', main = 'Below ground', cex.main = 2, cex.lab = 2, ...)
    axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1), tck = 0.015, 
      cex.axis = 1.5)
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 0), tck = 0.015, 
      cex.axis = 1.5)
    rect(xleft = test_points$mean_g - 0.505 * rect_size, xright = test_points$mean_g + 0.505 * rect_size, ybottom = test_points$pop - 0.505 * rect_size, 
      ytop = test_points$pop + 0.505 * rect_size, border = NA, col = adjustcolor(breakdown_col_pal$plow[pred_plow], alpha.f = 0.7))
    points(x = mean_g_scaled , y = pop_scaled, pch = 19, cex = 1.5, col = breakdown_col_pal$plow[action_space[Q$sample_points_df$policy, 'plow']]) 
    #dens
    plot(x = 0:1, y = 0:1, type = 'n', xlim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size), ylim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size),
      bty = 'n', tck = 0.02, xaxt = 'n', yaxt = 'n', main = 'Fecundity control', cex.lab = 2, cex.main = 2, ...)
    axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1), tck = 0.015, 
      cex.axis = 1.5)
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 0), tck = 0.015, 
      cex.axis = 1.5)
    rect(xleft = test_points$mean_g - 0.505 * rect_size, xright = test_points$mean_g + 0.505 * rect_size, ybottom = test_points$pop - 0.505 * rect_size, 
      ytop = test_points$pop + 0.505 * rect_size, border = NA, col = adjustcolor(breakdown_col_pal$dens[pred_dens], alpha.f = 0.7))
    points(x = mean_g_scaled , y = pop_scaled, pch = 19, cex = 1.5, col = breakdown_col_pal$dens[action_space[Q$sample_points_df$policy, 'dens']]) 
    #mech
    plot(x = 0:1, y = 0:1, type = 'n', xlim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size), ylim = c(0 - 0.51 * rect_size, 1 + 0.51 * rect_size),
      bty = 'n', tck = 0.02, xaxt = 'n', yaxt = 'n', main = 'Mechanical control', cex.main = 2, cex.lab = 2, ...)
    axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1), tck = 0.015, 
      cex.axis = 1.5)
    axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(Q$sample_points_df$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 0), tck = 0.015, 
      cex.axis = 1.5)
    rect(xleft = test_points$mean_g - 0.505 * rect_size, xright = test_points$mean_g + 0.505 * rect_size, ybottom = test_points$pop - 0.505 * rect_size, 
      ytop = test_points$pop + 0.505 * rect_size, border = NA, col = adjustcolor(breakdown_col_pal$mech[pred_mech], alpha.f = 0.7))
    points(x = mean_g_scaled , y = pop_scaled, pch = 19, cex = 1.5, col = breakdown_col_pal$mech[action_space[Q$sample_points_df$policy, 'mech']]) 
 
    #legend
    plot(x = 0:1, y = 0:1, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    legend(x = 'center', legend = c('Crop choice:', 'wheat', 'alternative', 'fallow', 'Direct managment:', 'no action', 'take action'),
      fill = c(NA, breakdown_col_pal$crop, NA, breakdown_col_pal$herb), 
      cex = 2, bty = 'n', border = NA, text.font = c(2, rep(1, 3), 2, rep(1, 2)))
  dev.off()
  setwd(current_loc)
}

#plot the result of the policy in terms of value, mean_g and black grass populaiton for every sampled point 

#plot the trajectory of a population given a starting state and a policy, have a value surface as the backgropund, then have dot with the action 
#and an arrow to the point in the state space where that action takes you and a dot with action taken there. Also have a number of moves to simulate
#so I can build up paths in succesive plots, also pass a vector of starting points to allow the comparision of a couple of starting points on the same plot
plot_simulated_policy <-function(Q, sim_obj, sims_to_plot, table_row_height = 0.05, output_loc, output_name = 'policy_trace_plot.pdf', ...){
  layout_mat = rbind(c(0, 0.7, 0, 1), c(0.7, 1, 0, 1)) 
  #plot value surface to show the movment over
  value_pallet = heat_hcl(n = 300)
  upper_mean = max(Q$sample_points_df$mean_g)
  lower_mean = min(Q$sample_points_df$mean_g)
  upper_pop = max(Q$sample_points_df$pop)
  lower_pop = min(Q$sample_points_df$pop)
  
  eval_res = 200
  value_grid = data.frame(mean_g = rep.int(seq(lower_mean, upper_mean, length.out = eval_res), eval_res), pop = rep(seq(lower_pop, upper_pop, length.out = eval_res), each = eval_res))#make every combination of evaluation points
  
  int_loc = getwd()
  setwd(output_loc) 
  pdf(file = output_name, width = 16, height = 10.1)
    for(i in sims_to_plot){
      split.screen(layout_mat)
      screen(1)
	#blank plot to add things to
	plot(0:1, 0:1, bty = 'n', xaxt = 'n', yaxt = 'n', type = 'n', ylab = 'population', xlab = 'resistance', ...)
	#add value surface
	if(class(Q$value_surface) == 'gbm'){
	  image(matrix(predict(Q$value_surface, newdata = value_grid, n.trees = Q$value_surface$n.trees), ncol = eval_res, byrow = FALSE), xaxt = 'n',  
	    yaxt = 'n', col = value_pallet, add = TRUE)
	}else{
	  image(matrix(predict(Q$value_surface, value_grid), ncol = eval_res, byrow = FALSE), xaxt = 'n', xlab = 'resistance', yaxt = 'n', ylab = 'population', add = TRUE)
	}
	axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(value_grid$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1))
	axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(value_grid$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 1))
	
	#extrac the points for the simulated object 
	mean_g_seq = as.numeric(c(sim_obj[[i]]$start_point['mean_g'], sapply(sim_obj[[i]]$state_over_time, FUN = function(x) x$next_state$mean_top)))
	pop_seq = as.numeric(c(sim_obj[[i]]$start_point['pop'], sapply(sim_obj[[i]]$state_over_time, FUN = function(x) x$next_state$pop_top)))
	reward_seq = round(as.numeric(sapply(sim_obj[[i]]$state_over_time, FUN = function(x) x$reward)), 1)
	
	#scaled state space parameters to fit on plot which is scaled from 0-1 by defualt
	mean_g_seq = (mean_g_seq - lower_mean) / (upper_mean - lower_mean)
	pop_seq = (pop_seq - lower_pop) / (upper_pop - lower_pop)
  
	#plot the movement as points on the value surface 
	for(j in 2:length(mean_g_seq)) arrowLine(x0 = mean_g_seq[j - 1], y0 = pop_seq[j - 1], x1 = mean_g_seq[j], y1 = pop_seq[j], lwd = 2, arr.width = 0.4, arr.length = 0.5) #make the lines with arrows
	points(x = mean_g_seq, y = pop_seq, pch = 19, cex = 5, col = c(rep(grey(0.9), length(mean_g_seq) - 1), grey(0.5)))
	text(x = mean_g_seq[1:(length(mean_g_seq) - 1)], y = pop_seq[1:(length(pop_seq) - 1)], labels = reward_seq, cex = 1)
	
      screen(2)
	#put the action sequence 
	#find all the non-na's
	non_na_acts = sim_obj[[i]]$action_seq[!is.na(sim_obj[[i]]$action_seq[ ,'act_num']), 'act_num']
	#set up the plotting space
	par(mar = c(0.5, 0.5, 3.5, 0.5))
	plot(0:1, 0:1, type = 'n', bty = 'n', main = '', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
	rect(xleft = 0, xright = 1, ybottom = 1 - (table_row_height * (2 + length(non_na_acts))), ytop = 1, col = grey(0.95), border = NA)
	
	col_xs = seq(0, 1 - (1/6), length = 6)
	row_ys = seq(1, 1 - (table_row_height * (dim(sim_obj[[i]]$action_seq)[1] + 1)), length.out = dim(sim_obj[[i]]$action_seq)[1] + 1)
	half_col = 0.5 * (1 / 6)
	half_row = 0.5 * table_row_height
	
	#headings
	text(x = col_xs + half_col, y = row_ys[1] - half_row, labels = c('time', 'crop', 'herb.', 'mech.', 'plow', 'dens.'), cex = 1.5)
	#times
	text(x = col_xs[1] + half_col, y = row_ys[2:(length(non_na_acts) + 1)] - half_row, labels = seq_along(non_na_acts), cex = 1.2)
	#crop
	crop_opts = c('W', 'A', 'F')
	text(x = col_xs[2] + half_col, y = row_ys[2:(length(non_na_acts) + 1)] - half_row, labels = crop_opts[sim_obj[[i]]$action_seq[seq_along(non_na_acts), 'crop']])
	
	#use points to show if the rest were used or not
	act_cols = c(grey(0.95), grey(0))
	points(rep(col_xs[3] + half_col, length(non_na_acts)), y = row_ys[2:(length(non_na_acts) + 1)] - half_row, pch = 19, cex = 1.5, 
	  col = act_cols[sim_obj[[i]]$action_seq[seq_along(non_na_acts), 'herb']])
	
	points(rep(col_xs[4] + half_col, length(non_na_acts)), y = row_ys[2:(length(non_na_acts) + 1)] - half_row, pch = 19, cex = 1.5, 
	  col = act_cols[sim_obj[[i]]$action_seq[seq_along(non_na_acts), 'mech']])
	
	points(rep(col_xs[5] + half_col, length(non_na_acts)), y = row_ys[2:(length(non_na_acts) + 1)] - half_row, pch = 19, cex = 1.5, 
	  col = act_cols[sim_obj[[i]]$action_seq[seq_along(non_na_acts), 'plow']])
	
	points(rep(col_xs[6] + half_col, length(non_na_acts)), y = row_ys[2:(length(non_na_acts) + 1)] - half_row, pch = 19, cex = 1.5, 
	  col = act_cols[sim_obj[[i]]$action_seq[seq_along(non_na_acts), 'dens']])
	#put some white lines so earier to follow the rows of the tabels
	for(j in 2:length(row_ys)) lines(y = c(row_ys[j], row_ys[j]), x = c(0, 1), col = 'white')
      close.screen(all.screens = TRUE)
    }
  dev.off()
  setwd(int_loc)
}

#function to draw arrows how I like them
arrowLine <- function(x0, y0, x1, y1, arr.width = 0.5, arr.length = 0.3, ...){
  lines(c(x0, x1), c(y0, y1), ...)
  Ax = seq(x0, x1, length = 3)
  Ay = seq(y0, y1, length = 3)
  Arrows(x0 = x0, y0 = y0, x1 = Ax[2], y1 = Ay[2], arr.width = arr.width, arr.length = arr.length)
}

#make a simple pretty plot of the value surface
plot_value_surface <- function(Q, output_loc, output_name = 'value_surface.pdf', ...){
  layout_mat = rbind(c(0, 0.85, 0, 1), c(0.85, 1, 0, 1)) 
  #plot value surface to show the movment over
  value_pallet = heat_hcl(n = 360)
  upper_mean = max(Q$sample_points_df$mean_g)
  lower_mean = min(Q$sample_points_df$mean_g)
  upper_pop = max(Q$sample_points_df$pop)
  lower_pop = min(Q$sample_points_df$pop)
  
  eval_res = 300
  value_grid = data.frame(mean_g = rep.int(seq(lower_mean, upper_mean, length.out = eval_res), eval_res), pop = rep(seq(lower_pop, upper_pop, length.out = eval_res), each = eval_res))#make every combination of evaluation points
  
  int_loc = getwd()
  setwd(output_loc) 
  pdf(file = output_name, width = 12, height = 10)
    split.screen(layout_mat)
    screen(1)
      #blank plot to add things to
      plot(0:1, 0:1, bty = 'n', xaxt = 'n', yaxt = 'n', type = 'n', ...)
      #add value surface
      if(class(Q$value_surface) == 'gbm'){
	predicted_values = matrix(predict(Q$value_surface, newdata = value_grid, n.trees = Q$value_surface$n.trees), ncol = eval_res, byrow = FALSE)
	image(predicted_values, xaxt = 'n', yaxt = 'n', col = value_pallet, add = TRUE)
      }else{
	predicted_values = matrix(predict(Q$value_surface, value_grid), ncol = eval_res, byrow = FALSE)
	image(predicted_values, xaxt = 'n', yaxt = 'n', col = value_pallet, add = TRUE)
      }
      axis(side = 1, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(value_grid$mean_g, probs = c(0, 0.25, 0.5, 0.75, 1)), 1), cex.axis = 1.5, tck = 0.015)
      axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(value_grid$pop, probs = c(0, 0.25, 0.5, 0.75, 1)), 0), cex.axis = 1.5, tck = 0.015)
      
    screen(2)
      #make a legend to show the scale 
      dummy_values = seq(min(predicted_values), max(predicted_values), length.out = 300)
      par(mar = c(3, 3, 3, 1))
      plot(0:1, 0:1, bty = 'n', xaxt = 'n', yaxt = 'n', type = 'n', xlab = '', ylab = '', main = 'value')
      image(t(as.matrix(dummy_values)), xaxt = 'n', yaxt = 'n', col = value_pallet, add = TRUE)
      axis(side = 2, at = c(0, 0.25, 0.5, 0.75, 1), labels = round(quantile(dummy_values, probs = c(0, 0.25, 0.5, 0.75, 1)), 0), las = 1)
    close.screen(all.screens = TRUE)
  dev.off()
  setwd(int_loc)
}
