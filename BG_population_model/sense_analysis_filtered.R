# sensitivity analysis of the 1D met and TS resistance model using meta-modeling
library(gbm)
library(gridExtra)

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

TODO: put the correlation plots of model summaries along with matrix plot of predictors, and BRT sense output for each measure
in a PDF, also fit the BRTs for each measure, so should be a 5 page summary pdf.

# check for relationship between the predictors

plot(sen_dat[, c('int_Rr', 'germ_prob', 'fec0', 'fec_cost', 'fec_max', 'dd_fec', 'herb_effect', 'g_prot', 'seed_sur', 
  'pro_exposed', 'scale_pollen', 'shape_pollen', 'seed_pro_short', 'seed_mean_dist_short', 'pro_seeds_to_mean_short', 
  'seed_mean_dist_long', 'pro_seeds_to_mean_long', 'final_ab_pop', 'final_x_occ', 'R_50', 'mean_spread', 'pro_rr_sur')])

#fit a BRT to explore which predictros affect different aspects of model behaviour  
BRT_pro_rr = gbm(qlogis(pro_rr_sur) ~ int_Rr + germ_prob + fec0 + fec_cost + fec_max + dd_fec + herb_effect + rel_g_pro + seed_sur +
  pro_exposed + scale_pollen + shape_pollen + seed_pro_short + seed_mean_dist_short + pro_seeds_to_mean_short + seed_mean_dist_long + pro_seeds_to_mean_long, 
  distribution = 'gaussian', interaction.depth = 4, shrinkage = 0.05, n.trees = 50000, cv.folds = 6, class.stratify.cv = TRUE, data = sen_dat, n.cores = 3, verbose = TRUE) 

# error with shrinkage = 0.05, op_trees = 20000 
# BRT_pro_rr$cv.error[op_trees]
# 0.04310812

  
  
# setwd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output")
# save(BRT_pro_rr, file = 'BRT_sense_pro_rr.Rdata') 
# load('BRT_sense_pro_rr.Rdata')

# extract useful info from the trees 
op_trees = gbm.perf(BRT_pro_rr, oobag.curve = TRUE, method = 'cv')

# get realtive influence
pro_rr_sum = summary(BRT_pro_rr, n.trees = op_trees, plotit = FALSE)

# looks like factors driving the differeence in pro_rr
# fec_cost, g_prot, int_Rr and pro_exposed (indexes: 4, 8, 1, 10)

plot_inds = c(4, 8, 1, 10)
plot_list = list()
par_sym = c('int[Rr]', 'phi[e]', 'f[0]', 'f[r]', 'f[max]', 'f[d]', 'xi','rho / xi', 'phi[b]', 'varsigma', 'a', 'c', 'alpha', 'mu[1]',
  'omega[1]', 'mu[2]', 'omega[2]')
#make a table to print 
rel_inf_df = data.frame(parameter = c( ,
   'phi[e]', 'a', seed_pro_short, shape_pollen, seed_mean_dist_short,  rep('a', 8)), 
  rel_inf = pro_rr_sum$rel.inf)
plot(tableGrob(rel_inf_df, cols = c("parameter", "rel_inf"), theme = ttheme_default(base_size = 9, parse = TRUE)))

plot_list[[1]] = tableGrob(rel_inf_df, cols = c("parameter", "rel_inf"), theme = ttheme_default(base_size = 9, parse = TRUE))
count = 2
for(i in 1:(length(plot_inds) - 1)){
  for(j in (i + 1):length(plot_inds)){
    plot_list[[count]] = plot(BRT_pro_rr, i.var = c(plot_inds[i], plot_inds[j]), type = 'response', n.trees = op_trees)
    count = count + 1
  }
}
#

setwd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output")
pdf(file = "sense_pro_rr.pdf", width = 10, height = 14.4)
  grid.arrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]], 
    plot_list[[6]], plot_list[[7]], ncol = 2)
dev.off()
