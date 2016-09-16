# sanity check script that looks for parameter combinations that produce population numbers 
# that are too low or too high compaitred to field observations. 
#install.packages('gbm', repos = 'http://cran.us.r-project.org')
library(gbm)
library(gridExtra)
library(gtable)
library(grid)
# get the data produced by the model runs
all_dat = read.csv("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/param_filtering_out.csv",
  header = TRUE, stringsAsFactors = FALSE)

# use post-herb above ground population to see select the parameter values that make 
# result in populations in the range 16,348 - 132,003, which comes from 
# Queenborough et al 2011, Figure 3 after some proccessing to get from the counts used in
# that study to plants per hectare. 
all_dat$final_pop = NA
all_dat$final_pop[all_dat$num_ab_ph_tot < 16348] = 'low'
all_dat$final_pop[all_dat$num_ab_ph_tot > 132000] = 'high'
all_dat$final_pop[is.na(all_dat$final_pop)] = 'in'
all_dat$in_out = ifelse(all_dat$final_pop == 'in', 1, 0)

# histograms of parameters for in vs out parameter values
par(mfrow = c(6, 3))
pred_names = names(all_dat)[1:30]
plot_inds = c(13, 15:30)
for(i in plot_inds){
  hist(all_dat[, i], main = pred_names[i], breaks = seq(min(all_dat[,i]), max(all_dat[,i]), length = 50))
  hist(all_dat[all_dat$in_out == 1, i], col = grey(0.5), border = grey(0.50), add = TRUE, 
    breaks = seq(min(all_dat[,i]), max(all_dat[,i]), length = 50))
}

# looks like it is almost all fec_max and fec_dd that controls if a population goes through the 
# sanity check or not
#Use a BRT to look for more complicated relationships that need intgeractions to explain them.
BRT_bi = gbm(in_out ~ int_Rr + germ_prob + fec0 + fec_cost + fec_max + dd_fec + herb_effect + g_prot + seed_sur +
  pro_exposed + scale_pollen + shape_pollen + seed_pro_short + seed_mean_dist_short + pro_seeds_to_mean_short + 
  seed_mean_dist_long + pro_seeds_to_mean_long, distribution = 'bernoulli', interaction.depth = 4, shrinkage = 0.05,
  n.trees = 10000, cv.folds = 10, class.stratify.cv = TRUE, data = all_dat, n.cores = 3, verbose = TRUE) 

# setwd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output")
# save(BRT_bi, file = 'BRT_pop_filter.Rdata') 
# load('BRT_pop_filter.Rdata')

# extract useful info from the trees 
op_trees = gbm.perf(BRT_bi, oobag.curve = TRUE, method = 'cv')

# get realtive influence
rel_inf = summary(BRT_bi, n.trees = op_trees)

#                                      var    rel.inf
# fec_max                           fec_max 41.9977878
# dd_fec                             dd_fec 34.3847773
# seed_sur                         seed_sur  7.2111341
# fec_cost                         fec_cost  4.7544122
# fec0                                 fec0  1.5901951
# germ_prob                       germ_prob  1.5785734
# seed_pro_short             seed_pro_short  1.2707816
# herb_effect                   herb_effect  1.1272816
# int_Rr                             int_Rr  1.1073332
# g_prot                             g_prot  1.0508439
# pro_exposed                   pro_exposed  0.9997911
# shape_pollen                 shape_pollen  0.9864782
# scale_pollen                 scale_pollen  0.9747225
# seed_mean_dist_short seed_mean_dist_short  0.9658881
# 

# looks like 2 important variables fec_max, and dd_fec, with 2 others varables, seed_sur and fec_cost,
# have some influence 
plot_inds = c(5, 6, 9, 4)
plot_list = list()
count = 1
for(i in 1:(length(plot_inds) - 1)){
  for(j in (i + 1):length(plot_inds)){
    plot_list[[count]] = plot(BRT_bi, i.var = c(plot_inds[i], plot_inds[j]), type = 'response', n.trees = op_trees)
    count = count + 1
  }
}
# looks like the fec_dd term could be a bit bigger maybe up to 0.15
setwd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output")
pdf(file = 'sanity_check_PDP.pdf', width = 10, height = 15)
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], ncol = 2)
dev.off()

passed_dat = all_dat[all_dat$in_out == 1, ]
write.csv(passed_dat, file = 'sanity_check_pars_passed.csv')


#make a table and plots of the relative influence and PDP for the 4 most important parameters in predicting in or out of the sanity check
par_sym_str = c('int[Rr]', 'phi[e]', 'f[0]', 'f[r]', 'f[max]', 'f[d]', 'xi', 'rho', 'phi[b]', 'varsigma', 'a', 'c', 'alpha', 'mu[1]',
  'omega[1]', 'mu[2]', 'omega[2]')
  
  
par_sym = c(expression(int[Rr]), expression(phi[e]), expression(f[0]), expression(f[r]), expression(f[max]), expression(f[d]), 
  expression(xi), expression(rho), expression(phi[b]), expression(varsigma), 'a', 'c', expression(alpha), expression(mu[1]),
  expression(omega[1]), expression(mu[2]), expression(omega[2]))

par_names = strsplit('int_Rr + germ_prob + fec0 + fec_cost + fec_max + dd_fec + herb_effect + g_prot + seed_sur + pro_exposed + scale_pollen + shape_pollen + seed_pro_short + seed_mean_dist_short + pro_seeds_to_mean_short + seed_mean_dist_long + pro_seeds_to_mean_long',
  split = ' + ', fixed = TRUE)[[1]]  
padding <- unit(5,"mm")
  
par_order = as.numeric(sapply(as.character(rel_inf$var), FUN = function(x) which(par_names == x)))
rel_inf_df = data.frame(parameter = par_sym_str[par_order], rel_inf = rel_inf$rel.inf)
table = tableGrob(rel_inf_df, cols = c("parameter", "rel. inf."), theme = ttheme_default(base_size = 6, parse = TRUE))

pdf(file = 'sanity_check_rel_inf.pdf', width = 2, height = 4)
  grid.draw(table)
dev.off()

plot_list = list()
plot_inds = c(5, 6, 9, 4)
count = 1
for(i in 1:(length(plot_inds) - 1)){
  for(j in (i + 1):length(plot_inds)){
    pg =  plot(BRT_bi, i.var = c(plot_inds[i], plot_inds[j]), n.trees = op_trees, return.grid = TRUE, type = 'response')
    preds = names(pg)
    form = paste0(preds[3], '~', preds[1], '+', preds[2])
    plot_list[[count]] = levelplot(as.formula(form), data = pg, xlab = list(label = par_sym[plot_inds[i]], cex = 1.7), 
      ylab = list(label = par_sym[plot_inds[j]], cex = 1.7), scales = list(x = list(cex = 1.5), y = list(cex = 1.5)), 
      colorkey = list(labels = list(cex = 1.5)))
    
    count = count + 1
  }
}

setwd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output")

pdf(file = 'sanity_check_PDP.pdf', width = 12, height = 15)
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], ncol = 2)
  grid.text(label = paste0(letters[1:6], ')'), x = c(0.05, 0.55), y = c(0.99, 0.99, 0.666, 0.666, 0.333, 0.333), gp = gpar(fontsize = 20))
dev.off()
