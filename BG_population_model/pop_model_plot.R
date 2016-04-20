# Functions to visulise the population model, filtering and sensitivity analysis.
library(colorspace)

## GET_MEAN_SD()
## takes a distribution and gets the mean, sd and total number of seeds in the seedbank
get_mean_sd <- function(dist, eval_points, dg){
  total_sum = sum(dist) * dg
  approx_mean = sum((dist / total_sum) * eval_points * dg)
  approx_sd = sqrt(sum(((eval_points - approx_mean)^2) * (dist / total_sum) * dg))
  return(list(approx_mean = approx_mean, approx_sd = approx_sd, total_pop = total_sum))
}

## POP_PLOTTER()
## takes a list of populatin matricies over time and 
pop_plotter <- function(pop_list, eval_points, dg, color_pallet = rainbow_hcl(length(pop_list)), 
  line_labels = paste0('line', 1:length(pop_list)), out_name = 'generic_plot.pdf', out_loc){
  
  n_time = dim(pop_list[[1]])[1]
  pop = matrix(NA, nrow = length(pop_list), ncol = n_time)
  mean_g = matrix(NA, nrow = length(pop_list), ncol = n_time)
  sd_g = matrix(NA, nrow = length(pop_list), ncol = n_time)
  
  #get the populations for each population in the list
  for(i in 1:length(pop_list)){
    ncols = dim(pop_list[[i]])[2]
    out_list = apply(pop_list[[i]], MARGIN = 1, FUN = function(x){
      get_mean_sd(x, eval_points, dg)
    })
    pop[i, ] = sapply(out_list, FUN = function(x) x$total_pop)
    mean_g[i, ] = sapply(out_list, FUN = function(x) x$approx_mean)
    sd_g[i, ] = sapply(out_list, FUN = function(x) x$approx_sd)
  }
  
  int_loc = getwd()
  setwd(out_loc)
  pdf(file = out_name, width = 10, height = 10)
    par(mfrow = c(2, 2))
    #plot of population over time
    plot(1:n_time, pop[1, ], type = 'n', bty = 'n', xlab = 'time', ylab = 'population', 
      ylim = c(0, max(pop)))
    for(i in 1:dim(pop)[1]) lines(1:n_time, pop[i, ], lwd = 1.5, col = color_pallet[i]) 
    legend(x = 'bottomright', legend = line_labels, fill = color_pallet, border = NA, bty = 'n')
    
    #plot of mean g over time
    plot(1:n_time, mean_g[1, ], type = 'n', bty = 'n', xlab = 'time', ylab = 'mean_g', 
      ylim = c(min(mean_g), max(mean_g)))
    for(i in 1:dim(mean_g)[1]) lines(1:n_time, mean_g[i, ], lwd = 1.5, col = color_pallet[i]) 
    
    #plot of sd over time
    plot(1:n_time, sd_g[1, ], type = 'n', bty = 'n', xlab = 'time', ylab = 'sd_g', 
      ylim = c(min(sd_g), max(sd_g)))
    for(i in 1:dim(sd_g)[1]) lines(1:n_time, sd_g[i, ], lwd = 1.5, col = color_pallet[i]) 
    
    #plot of intial populaiton, mid-time populaiton and end population over g
    temp_pallet = coords(as(hex2RGB(color_pallet), 'polarLUV'))
    temp_pallet[, "C"] <- 25
    color_pallet2 = hex(polarLUV(temp_pallet), fixup = TRUE)
    
    max_pop_g = max(sapply(pop_list, FUN = max))
    
    plot(x = eval_points, y = seq(0, max_pop_g, length = length(eval_points)), type = 'n', 
      bty = 'n', xlab = 'resistance (g)', ylab = 'population', log = 'y')
    for(i in 1:length(pop_list)){
      #first plot the final pop distribution
      polygon(x = eval_points, y = pop_list[[i]][dim(pop_list[[i]])[2], ], col = color_pallet[i], 
	border = NA)
      #then plot the first population
      polygon(x = eval_points, y = pop_list[[i]][1, ], col = color_pallet2[i], border = NA)
    }
  dev.off()
  setwd(int_loc)
}


## TSR_POP_v_time_PLOT()
## takes the output of a pop_run_TSR() and make a plot of the populaiton changes over time

TSR_pop_v_time <- function(pops_list, plot_inds, out_loc, out_name = 'pop_v_time_TSR_all.pdf'){

  layout_mat = rbind(c(0, 1, 0.6, 1), c(0, 1, 0.4,  0.6), c(0, 1, 0, 0.4))
  
  int_loc = getwd()
  setwd(out_loc)
  pdf(file = 'filtered_pops_plots.pdf', height = 20, width = 10)
    for(i in plot_inds){
      split.screen(layout_mat)
      #above ground plot
      screen(1)
	max_y = max(c(pops_list$ag_pre_herb[, i], pops_list$ag_noherb[, i], pops_list$ag_post_herb[, i]))
	plot(x = 1:dim(pops_list$ag_noherb)[1], y = pops_list$ag_noherb[, i], type = 'l', col = 'blue', bty = 'n', 
	  ylim = c(0, max_y), ylab = 'above ground population', xlab = '', main = 'above ground', lwd = 2,
	  cex.main = 2, cex.lab = 2, cex.axis = 1.7)
	lines(x = 1:dim(pops_list$ag_noherb)[1], y = pops_list$ag_pre_herb[, i], col = 'green', lwd = 2)
	lines(x = 1:dim(pops_list$ag_noherb)[1], y = pops_list$ag_post_herb[, i], col = 'red', lwd = 2)
      
      #TSR plot
      screen(2)
	plot(x = 1:dim(pops_list$TSR_freq[[i]])[1], y = seq(0, 1, length = dim(pops_list$TSR_freq[[i]])[1]), 
	  type = 'n', ylim = c(0, 1), xlim = c(1, dim(pops_list$TSR_freq[[i]])[1]), bty = 'n', xlab = '',
	  ylab = 'mean_g', main = 'TS and metabolic resistance', yaxt = 'n', tck = 0.02, cex.lab = 2, 
	  cex.axis = 1.7, cex.main = 2)
	polygon(x = c(0.99999, 1:dim(pops_list$TSR_freq[[i]])[1], dim(pops_list$TSR_freq[[i]])[1] + 1.0e-05), 
	  y = c(0, rep(1, dim(pops_list$TSR_freq[[i]])[1]), 0), col = 'skyblue', border = NA)
	polygon(x = c(0.99999, 1:dim(pops_list$TSR_freq[[i]])[1], dim(pops_list$TSR_freq[[i]])[1] + 1.0e-05), 
	  y = c(0, rowSums(pops_list$TSR_freq[[i]][, 1:2]) , 0), col = 'magenta', border = NA)
	polygon(x = c(0.99999, 1:dim(pops_list$TSR_freq[[i]])[1], dim(pops_list$TSR_freq[[i]])[1] + 1.0e-05), 
	  y = c(0, pops_list$TSR_freq[[i]][, 1], 0), col = 'red', border = NA)
	#add the maen_g plot over the top, scalling the max to 0.9 and min to 0.1 
	range_g = range(pops_list$seedbank_g_herb[, i])
	rescaled_g = 0.1 + (0.9 - 0.1) * (pops_list$seedbank_g_herb[, i] - range_g[1]) /
	  (range_g[2] - range_g[1])
	lines(x = 1:dim(pops_list$TSR_freq[[i]])[1], y = rescaled_g, lwd = 2)
	axis(side = 2, at = c(0.1, 0.5, 0.9), labels = round(c(range_g[1], ((range_g[2] - range_g[1]) / 2) + range_g[1],
	  range_g[2]), 4), tck = 0.02, cex.axis = 1.7)	
	  
      #seedbank plot 
      screen(3)
	max_y = max(c(pops_list$seedbank_herb[, i], pops_list$seedbank_noherb[, i]))
	plot(x = 1:dim(pops_list$seedbank_herb)[1], y = pops_list$seedbank_herb[, i], type = 'l', col = 'red', bty = 'n', 
	  ylim = c(0, max_y), ylab = 'seedbank population', xlab = 'time', main = 'seedbank', lwd = 2, 
	  cex.main = 2, cex.lab = 2, cex.axis = 1.7, tck = 0.02)
	lines(x = 1:dim(pops_list$seedbank_herb)[1], y = pops_list$seedbank_noherb[, i], col = 'blue', lwd = 2)
      close.screen(all.screens = TRUE)
    }
  dev.off()
  setwd(int_loc)
}



