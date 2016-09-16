# sensitivity analysis of the 1D met and TS resistance model using meta-modeling
library(gbm)
library(gridExtra)
library(grid)
library(gtable)
library(colorspace)
library(hexbin)

sen_dat = read.csv("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/senes_runs.csv",
  header = TRUE, stringsAsFactors = FALSE)

sen_dat$rel_g_pro = sen_dat$g_prot / sen_dat$herb_effect
sen_dat$max_g_sur = 1 / (1 + exp(-(sen_dat$base_sur - (sen_dat$herb_effect - 
  pmin(sen_dat$herb_effect, sen_dat$g_prot * sen_dat$max_g)))))  
  
# check the relationship between the measures of pop_performance
plot(sen_dat[, c('final_ab_pop', 'final_x_occ', 'final_R', 'R_50', 'mean_spread', 'pro_rr_sur', 'time_2_max_g', 
  'max_g', 'final_g', 'max_g_sur')]) #add others from new runs

# final_x_occ is 1:1 correlated with mean_spread also R_50 and final_R are closley correlated since most runs did not reach 100% R. 
# final_R and pro_rr_sur have a very tight correlation, as one likely blocks the other so best 3 to look at are 
# final_R, mean_spread and max_g_sur
setwd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output")

BRT_fin_R = gbm(final_R ~ int_Rr + germ_prob + fec0 + fec_cost + fec_max + dd_fec + herb_effect + rel_g_pro + seed_sur +
  pro_exposed + scale_pollen + shape_pollen + seed_pro_short + seed_mean_dist_short + pro_seeds_to_mean_short + seed_mean_dist_long + pro_seeds_to_mean_long, 
  distribution = 'gaussian', interaction.depth = 4, shrinkage = 0.05, n.trees = 50000, cv.folds = 6, data = sen_dat, n.cores = 3, verbose = TRUE) 
save(BRT_fin_R, file = 'BRT_sense_fin_R.Rdata') 
# load('BRT_sense_fin_R.Rdata')

BRT_mean_spread = gbm(mean_spread ~ int_Rr + germ_prob + fec0 + fec_cost + fec_max + dd_fec + herb_effect + rel_g_pro + seed_sur +
  pro_exposed + scale_pollen + shape_pollen + seed_pro_short + seed_mean_dist_short + pro_seeds_to_mean_short + seed_mean_dist_long + pro_seeds_to_mean_long, 
  distribution = 'gaussian', interaction.depth = 4, shrinkage = 0.05, n.trees = 50000, cv.folds = 6, data = sen_dat, n.cores = 3, verbose = TRUE) 
save(BRT_mean_spread, file = 'BRT_sense_mean_spread.Rdata') 
# load('BRT_sense_mean_spread.Rdata')
  
BRT_max_g_sur = gbm(max_g_sur ~ int_Rr + germ_prob + fec0 + fec_cost + fec_max + dd_fec + herb_effect + rel_g_pro + seed_sur +
  pro_exposed + scale_pollen + shape_pollen + seed_pro_short + seed_mean_dist_short + pro_seeds_to_mean_short + seed_mean_dist_long + pro_seeds_to_mean_long, 
  distribution = 'gaussian', interaction.depth = 4, shrinkage = 0.05, n.trees = 50000, cv.folds = 6, data = sen_dat, n.cores = 3, verbose = TRUE) 
save(BRT_max_g_sur, file = 'BRT_sense_max_g_sur.Rdata') 
# load('BRT_sense_max_g_sur.Rdata')

op_trees_fR = gbm.perf(BRT_fin_R, oobag.curve = TRUE, method = 'cv')
op_trees_ms = gbm.perf(BRT_mean_spread, oobag.curve = TRUE, method = 'cv')
op_trees_mg = gbm.perf(BRT_max_g_sur, oobag.curve = TRUE, method = 'cv')

BRT_fin_R$cv.error[op_trees_fR] # 0.001196987
BRT_mean_spread$cv.error[op_trees_ms] # 1.912273e-08
BRT_mean_spread$cv.error[op_trees_mg] # 2.493968e-08

sum_fR = summary(BRT_fin_R, n.trees = op_trees_fR , plotit = FALSE)
sum_ms = summary(BRT_mean_spread, n.trees = op_trees_ms, plotit = FALSE)
sum_mg = summary(BRT_max_g_sur, n.trees = op_trees_mg, plotit = FALSE)

# make the tables of relative influence for each measure
par_sym = c(expression(int[Rr]), expression(phi[e]), expression(f[0]), expression(f[r]), expression(f[max]), expression(f[d]), 
  expression(xi), expression(rho / xi), expression(phi[b]), expression(varsigma), 'a', 'c', expression(alpha), expression(mu[1]),
  expression(omega[1]), expression(mu[2]), expression(omega[2]))
par_sym_str = c('int[Rr]', 'phi[e]', 'f[0]', 'f[r]', 'f[max]', 'f[d]', 'xi', 'rho / xi', 'phi[b]', 'varsigma', 'a', 'c', 'alpha', 'mu[1]',
  'omega[1]', 'mu[2]', 'omega[2]')
  
  
par_names = strsplit('int_Rr + germ_prob + fec0 + fec_cost + fec_max + dd_fec + herb_effect + rel_g_pro + seed_sur + pro_exposed + scale_pollen + shape_pollen + seed_pro_short + seed_mean_dist_short + pro_seeds_to_mean_short + seed_mean_dist_long + pro_seeds_to_mean_long',
  split = ' + ', fixed = TRUE)[[1]]  
padding <- unit(5,"mm")
  
fR_par_order = sapply(as.character(sum_fR$var), FUN = function(x) which(par_names == x))
rel_inf_df_fR	 = data.frame(parameter = par_sym_str[fR_par_order], 
  rel_inf = sum_fR$rel.inf)
title <- textGrob("final R", gp = gpar(fontsize = 15))
table_0 = tableGrob(rel_inf_df_fR, cols = c("parameter", "rel_inf"), theme = ttheme_default(base_size = 9, parse = TRUE))
table = gtable_add_rows(table_0, heights = grobHeight(title) + padding, pos = 0)
fR_table = gtable_add_grob(table, title, 1, 1, 1, ncol(table_0))

ms_par_order = sapply(as.character(sum_ms$var), FUN = function(x) which(par_names == x))
rel_inf_df_ms	 = data.frame(parameter = par_sym_str[ms_par_order], 
  rel_inf = sum_ms$rel.inf)
title <- textGrob("mean spread", gp = gpar(fontsize = 15))
table_0 = tableGrob(rel_inf_df_ms, cols = c("parameter", "rel_inf"), theme = ttheme_default(base_size = 9, parse = TRUE))
table = gtable_add_rows(table_0, heights = grobHeight(title) + padding, pos = 0)
ms_table = gtable_add_grob(table, title, 1, 1, 1, ncol(table_0))

mg_par_order = sapply(as.character(sum_mg$var), FUN = function(x) which(par_names == x))
rel_inf_df_mg	 = data.frame(parameter = par_sym_str[mg_par_order], 
  rel_inf = sum_mg$rel.inf)
title <- textGrob("sur max g", gp = gpar(fontsize = 15))
table_0 = tableGrob(rel_inf_df_mg, cols = c("parameter", "rel_inf"), theme = ttheme_default(base_size = 9, parse = TRUE))
table = gtable_add_rows(table_0, heights = grobHeight(title) + padding, pos = 0)
mg_table = gtable_add_grob(table, title, 1, 1, 1, ncol(table_0))

# plot the tables
pdf('sense_rel_inf_tables.pdf', height = 5.5, width = 6)
  grid.arrange(fR_table, ms_table, mg_table, ncol = 3, widths = c(1, 1, 1), heights = c(0.2))
dev.off()

 #PDP plots for final_R
pdf(file = 'sense_PDP_fin_R.pdf', height = 15, width = 12)
  plot_inds = c(4, 8, 1, 3)
  plot_list = list()
  count = 1
  for(i in 1:(length(plot_inds) - 1)){
    for(j in (i + 1):length(plot_inds)){
      pg =  plot(BRT_fin_R, i.var = c(plot_inds[i], plot_inds[j]), n.trees = op_trees_fR, return.grid = TRUE)
      preds = names(pg)
      form = paste0(preds[3], '~', preds[1], '+', preds[2])
      plot_list[[count]] = levelplot(as.formula(form), data = pg, xlab = list(label = par_sym[plot_inds[i]], cex = 1.7), 
	ylab = list(label = par_sym[plot_inds[j]], cex = 1.7), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), 
	colorkey = list(labels = list(cex = 1.5)))
      
      count = count + 1
    }
  }
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], 
    plot_list[[6]], ncol = 2)#, top = textGrob('final R', gp = gpar(cex = 3)))
  grid.text(label = paste0(letters[1:6], ')'), x = c(0.05, 0.55), y = c(0.99, 0.99, 0.666, 0.666, 0.333, 0.333), gp = gpar(fontsize = 20))
dev.off()

#PDP for mean spread, top 5 predictors
pdf(file = 'sense_PDP_mean_spread.pdf', height = 15, width = 15)
  plot_inds = c(17, 5, 9, 16, 6)
  plot_list = list()
  count = 1
  for(i in 1:(length(plot_inds) - 1)){
    for(j in (i + 1):length(plot_inds)){
      pg =  plot(BRT_mean_spread, i.var = c(plot_inds[i], plot_inds[j]), n.trees = op_trees_fR, return.grid = TRUE)
      preds = names(pg)
      form = paste0(preds[3], '~', preds[1], '+', preds[2])
      plot_list[[count]] = levelplot(as.formula(form), data = pg, xlab = list(label = par_sym[plot_inds[i]], cex = 1.7), 
	ylab = list(label = par_sym[plot_inds[j]], cex = 1.7), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), 
	colorkey = list(labels = list(cex = 1.5)))
      
      count = count + 1
    }
  }
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], 
    plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], plot_list[[10]], 
    ncol = 3)#, top = textGrob('mean spread', gp = gpar(cex = 3)))
  grid.text(label = paste0(letters[1:10], ')'), x = c(0.05, 0.35, 0.68), y = c(0.99, 0.99, 0.99, 0.75, 0.75, 0.75, 0.5, 0.5, 0.5, 0.25), gp = gpar(fontsize = 20))
dev.off()

#PDP for max_sur_g
pdf(file = 'sense_PDP_max_g_sur.pdf', height = 15, width = 12)
  plot_inds = c(8, 4, 1, 10)
  plot_list = list()
  count = 1
  for(i in 1:(length(plot_inds) - 1)){
    for(j in (i + 1):length(plot_inds)){
      pg =  plot(BRT_max_g_sur, i.var = c(plot_inds[i], plot_inds[j]), n.trees = op_trees_fR, return.grid = TRUE)
      preds = names(pg)
      form = paste0(preds[3], '~', preds[1], '+', preds[2])
      plot_list[[count]] = levelplot(as.formula(form), data = pg, xlab = list(label = par_sym[plot_inds[i]], cex = 1.7), 
	ylab = list(label = par_sym[plot_inds[j]], cex = 1.7), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), 
	colorkey = list(labels = list(cex = 1.5)))
      
      count = count + 1
    }
  }
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], 
    plot_list[[6]], ncol = 2)#, top = textGrob('survival under max g', gp = gpar(cex = 3)))
  grid.text(label = paste0(letters[1:6], ')'), x = c(0.05, 0.55), y = c(0.99, 0.99, 0.666, 0.666, 0.333, 0.333), gp = gpar(fontsize = 20))
dev.off()

# plot to show the relationship between max g, time to max g and the final R, showing interaction between 
# the two pathways to resistance.
fin_R_col = sequential_hcl(n = 100, h = c(260, 360), c = 80, l = 90, alpha = 0.2)

fr_col_fun = colorRamp(c('red', 'blue'), space = 'rgb')
fr_col_raw = fr_col_fun(sen_dat$final_R)
fr_col_raw[, 2] = fr_col_raw[, 2] + 0.0001
fin_R_col = adjustcolor(fr_col_raw, alpha.f = 0.2)

plot(sen_dat$max_g_sur, sen_dat$time_2_max_g, bty = 'n', tck = 0.015, col = fin_R_col)
plot(sen_dat$final_R, sen_dat$final_R, bty = 'n', tck = 0.015, col = fin_R_col)


fr_col_fun = colorRamp(c('red', 'blue'), space = 'rgb', interpolate = 'linear')
fr_col_raw = fr_col_fun(seq(0, 1, 0.1))
fr_col_raw[, 2] = fr_col_raw[, 2] + 0.0001
fin_R_col = adjustcolor(fr_col_raw, alpha.f = 0.7)
plot(seq(0, 1, 0.1), col = fr_col_raw)


# better 3D hist function to let me use two color channels so intensity is the count and hue is a variable
# counts made on x and y (color intensity), hue codes z
hist_3D = function(x, y, z, nbins = 30, min_chrom = 5, start_hue = 240, end_hue = 360, ...){
  #bins for counts
  xbins = seq(min(x), max(x), length = nbins) 
  ybins = seq(min(y), max(y), length = nbins)
  
  bins = matrix(NA, ncol = 4, nrow = (length(xbins) - 1) * (length(ybins) - 1))
  count = 1 
  for(x_ind in 2:length(xbins)){
    for(y_ind in 2:length(ybins)){
      bins[count, ] = c(xbins[x_ind - 1], xbins[x_ind], ybins[y_ind - 1], ybins[y_ind])
      count = count + 1
    }
  }
  
  #get the count for each bin and the mean z
  bin_count = numeric(dim(bins)[1])
  bin_z = numeric(dim(bins)[1])
  bin_col = NA  
  for(i in 1:dim(bins)[1]){
    in_bin = x >= bins[i, 1] & x < bins[i, 2] & y >= bins[i, 3] & y < bins[i, 4]
    bin_count[i] = sum(in_bin)
    bin_z[i] = ifelse(bin_count[i] == 0, 0, mean(z[in_bin]))
  }
  
  bin_count = ifelse(bin_count == 0, 0, bin_count + min_chrom)
  z_min = min(bin_z[bin_z > 0])
  z_max = max(bin_z)
  hue = ifelse(bin_z == 0, 180, (((end_hue - start_hue) * (bin_z - z_min)) / (z_max - z_min)) + start_hue)
  c_min = min(bin_count[bin_count > 0])
  c_max = max(bin_count[bin_count > 0])
  chrom = ifelse(bin_count == 0, 0, ((((35 + min_chrom)  - min_chrom) * (bin_count - c_min)) / (c_max - c_min)) + min_chrom)
  bin_col = hcl(h = hue, c = chrom, l = 98)
  
  plot(seq(min(x), max(x), length = nbins), seq(min(y), max(y), length = nbins), type = 'n', bty = 'n', ...)
  rect(xleft = bins[, 1], ybottom = bins[, 3], xright = bins[, 2], ytop = bins[, 4], col = bin_col, border = bin_col)
}

# make a color legend for the above
hist_3D_leg = function(min_chrom = 5, start_hue = 240, end_hue = 360, ...){
  
  xbins = seq(0, 1, length = 50) 
  ybins = seq(0, 1, length = 50)
  
  hue_range = seq(0, 1, length = 49)
  co_range = seq(0, 1, length = 49)
  
  bins = matrix(NA, ncol = 6, nrow = (length(xbins) - 1) * (length(ybins) - 1))
  count = 1 
  for(x_ind in 2:length(xbins)){
    for(y_ind in 2:length(ybins)){
      bins[count, 1:4] = c(xbins[x_ind - 1], xbins[x_ind], ybins[y_ind - 1], ybins[y_ind])
      bins[count, 5] = hue_range[y_ind - 1]
      bins[count, 6] = co_range[x_ind - 1]
      count = count + 1
    }
  }
 
  colnames(bins) = c('xl', 'xr', 'yb', 'yt', 'hu', 'co')
  hue = ifelse(bins[, 'hu'] == 0, 180, (((end_hue - start_hue) * (bins[, 'hu'] - bins[2, 'hu'])) / (1 - bins[2, 'hu'])) + start_hue)
  chrom = ifelse(bins[, 'co'] == 0, 0, ((((35 + min_chrom)  - min_chrom) * (bins[, 'co'] - bins[50, 'co'])) / (1 - bins[50, 'co'])) + min_chrom)
  
  bin_col = hcl(h = hue, c = chrom, l = 98)
  plot(seq(0, 1, length = 50), seq(0, 1, length = 50), type = 'n', bty = 'n', ...)
  rect(xleft = bins[, 'xl'], ybottom = bins[, 'yb'], xright = bins[, 'xr'], ytop = bins[, 'yt'], col = bin_col, border = bin_col)
} 

  

hist_3D(x = sen_dat$max_g_sur, y = sen_dat$time_2_max_g, z = sen_dat$final_R, nbins = 15, min_chrom = 25, start_hue = 180, end_hue = 360, 
  xlab = 'survival under max g', ylab = 'time to max g', cex.lab = 1.5, tck = 0.015, cex.axis = 1.5)

hist_3D_leg(min_chrom = 25, start_hue = 180, end_hue = 360)

#does not show too much, just use a matrix of hex plots to show the relationships
plot(textGrob('final R', gp = gpar(cex = 3))
grid.text('final R')

finR_lab = grid.text("final R", x = 0.5, y = 0.5, rot = 0, gp = gpar(fontsize = 20), 
  check = FALSE, draw = FALSE)
mg_lab = grid.text("survival under\nmax g", x = 0.5, y = 0.5, rot = 0, gp = gpar(fontsize = 20), 
  check = FALSE, draw = FALSE)
t2mg_lab = grid.text("time to\nmax g", x = 0.5, y = 0.5, rot = 0, gp = gpar(fontsize = 20), 
  check = FALSE, draw = FALSE)
ms_lab = grid.text("mean\nspread", x = 0.5, y = 0.5, rot = 0, gp = gpar(fontsize = 20), 
  check = FALSE, draw = FALSE)
blank = grid.text("", x = 0.5, y = 0.5, rot = 0, gp = gpar(fontsize = 20), 
  check = FALSE, draw = FALSE)

  
mg_v_fR = hexbinplot(final_R ~ max_g_sur, data = sen_dat, xlab = list(label = ''), 
  ylab = list(label = ''), aspect = 1)
t2mg_v_fR = hexbinplot(final_R ~ time_2_max_g, data = sen_dat, xlab = list(label = ''), 
  ylab = list(label = ''), aspect = 1)
ms_v_fR = hexbinplot(final_R ~ mean_spread, data = sen_dat, xlab = list(label = ''), 
  ylab = list(label = ''), aspect = 1)

t2mg_v_mg = hexbinplot(max_g_sur ~ time_2_max_g, data = sen_dat, xlab = list(label = ''), 
  ylab = list(label = ''), aspect = 1)
ms_v_mg = hexbinplot(max_g_sur ~ mean_spread, data = sen_dat, xlab = list(label = ''), 
  ylab = list(label = ''), aspect = 1)

ms_v_t2mg = hexbinplot(time_2_max_g ~ mean_spread, data = sen_dat, xlab = list(label = ''), 
  ylab = list(label = ''), aspect = 1)

pdf('model_behave_plots.pdf', height = 20, width = 20)
  grid.arrange(finR_lab, mg_v_fR, t2mg_v_fR, ms_v_fR, blank, 
    mg_lab, t2mg_v_mg, ms_v_mg, blank, blank, 
    t2mg_lab, ms_v_t2mg, blank, blank, blank, ms_lab, 
    ncol = 4, widths = c(0.5, 1, 1, 1), 
    heights = c(1, 1, 1, 0.5)) 
dev.off()  





