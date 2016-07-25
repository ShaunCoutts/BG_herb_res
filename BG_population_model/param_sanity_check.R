# sanity check script that looks for parameter combinations that produce population numbers 
# that are too low or too high compaitred to field observations. 
#install.packages('gbm', repos = 'http://cran.us.r-project.org')
library(gbm)
library(gridExtra)

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
  hist(all_dat[, i], main = pred_names[i])
  hist(all_dat[all_dat$in_out == 'in', i], col = grey(0.5), border = grey(0.50), add = TRUE)
}

# looks like it is almost all fec_max and fec_dd that controls if a population goes through the 
# sanity check or not
#Use a BRT to look for more complicated relationships that need intgeractions to explain them.
BRT_bi = gbm(in_out ~ int_Rr + germ_prob + fec0 + fec_cost + fec_max + dd_fec + herb_effect + g_prot + seed_sur +
  pro_exposed + scale_pollen + shape_pollen + seed_pro_short + seed_mean_dist_short, + pro_seeds_to_mean_short + 
  seed_mean_dist_long + pro_seeds_to_mean_long, distribution = 'bernoulli', interaction.depth = 4, shrinkage = 0.05,
  n.trees = 10000, cv.folds = 10, class.stratify.cv = TRUE, data = all_dat, n.cores = 3, verbose = TRUE) 

# setwd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output")
# save(BRT_bi, file = 'BRT_pop_filter.Rdata') 
# load('BRT_pop_filter.Rdata')

# extract useful info from the trees 
op_trees = gbm.perf(BRT_bi, oobag.curve = TRUE, method = 'cv')

# get realtive influence
summary(BRT_bi, n.trees = op_trees)
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

grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], plot_list[[6]], ncol = 2)

# looks like the fec_dd term could be a bit bigger maybe up to 0.15

