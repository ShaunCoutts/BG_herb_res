# natural spread experiments 
#######################################################################################################
# make a wrapper to run the natural spread scenarios in parallel  

########################################################################################################3
using DataFrames
output_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/" 
@everywhere file_loc_func_p = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model" 
@everywhere cd(file_loc_func_p)
@everywhere include("BG_met_TSR_space_exper_funs.jl")
@everywhere include("BG_met_TSR_space_pop_process.jl")
@everywhere include("BG_met_TSR_space_dispersal_functions.jl")
include("spatial_model_plotting_script.jl")
  

###################################################################################################
# script to run the translocation experiments
# parameter values
upper_g = 10;
lower_g = -10;
dg = 0.5;
int_mean_g = 0.0;
int_sd_g = 1.4142;

x_dim = 150; # number of spatial evaluation points, actual landscape size is x_dim * dx
dx = 1.0;
source_locs = 1:20 # source population locations, at edge of landscape
recive_locs = 21:x_dim

int_num_RR = 0.1;
int_num_Rr = 0.0;
int_num_rr = 10.0; # number of intial seeds at each location for each genoptype, assume only TS susceptible
num_iter = 150;

seed_pro_short = 0.48; 
seed_mean_dist_short = 0.58; 
pro_seeds_to_mean_short = 0.44; 
seed_mean_dist_long = 1.65; 
pro_seeds_to_mean_long = 0.39;
scale_pollen = 32.0;
shape_pollen = 3.32; 
offspring_sd = 1.0;
fec_max = 60.0;
fec0 = 4.0;
fec_cost = 0.5;  # when fec0 = 5 and fec_cost = 0.8 fecundity is fec_max*0.62 at g = 5
dd_fec = 0.15;
base_sur = 10.0; 
herb_effect = 14.0; # herbicide kills 98% of suceptable plants
g_prot = 1.5; # when g = 5 survival will be 0.97 * base_sur (basically same as TSR). 
pro_exposed = 0.8;
seed_sur = 0.45;
germ_prob = 0.52;
resist_G = ["RR", "Rr"];

# set up the evaluation points for quantitative resistance
g_vals = collect(lower_g : dg : upper_g);   

# seed and pollen dispersal mixing kernels
seed_disp_mat_1D = zeros(x_dim, x_dim);
seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, seed_mean_dist_short, 
  pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long);

pollen_disp_mat = zeros(x_dim, x_dim);
pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx, a = scale_pollen, c = shape_pollen);

# effect of g on survival and fecundity
g_effect_fec = exp(-(fec0 - abs(g_vals) * fec_cost));
sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed);

# set up the different source/reciving configurations
source_rr_reciv_empty = zeros(x_dim);
source_rr_reciv_empty[source_locs] = int_num_rr; 
source_Rr_reciv_empty = zeros(x_dim);
source_Rr_reciv_empty[source_locs] = int_num_Rr; 
source_RR_reciv_empty = zeros(x_dim);
source_RR_reciv_empty[source_locs] = int_num_RR; 

source_rr_reciv_full = zeros(x_dim) + int_num_rr;
source_Rr_reciv_full = zeros(x_dim);
source_Rr_reciv_full[source_locs] = int_num_Rr; 
source_RR_reciv_full = zeros(x_dim);
source_RR_reciv_full[source_locs] = int_num_RR; 

# set up the herbicide application scenarios
herb_app_source_expos = convert(Array{Int64, 1}, ones(x_dim));
herb_app_source_expos[source_locs] += 1;

herb_app_reciv_expos = convert(Array{Int64, 1}, ones(x_dim));
herb_app_reciv_expos[recive_locs] += 1;

herb_app_all_expos = convert(Array{Int64, 1}, ones(x_dim)) + 1;
herb_app_none_expos = convert(Array{Int64, 1}, ones(x_dim));

# run all 8 scenarios to produce the time series for a plot
# herb_nn_ls_empty = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_empty, 
#   source_Rr_reciv_empty, source_RR_reciv_empty, int_mean_g, int_sd_g, seed_sur, 
#   germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
#   seed_disp_mat_1D, pollen_disp_mat, herb_app_none_expos);
# 
herb_en_ls_empty = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_empty, 
  source_Rr_reciv_empty, source_RR_reciv_empty, int_mean_g, int_sd_g, seed_sur, 
  germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
  seed_disp_mat_1D, pollen_disp_mat, herb_app_source_expos);

herb_ne_ls_empty = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_empty, 
  source_Rr_reciv_empty, source_RR_reciv_empty, int_mean_g, int_sd_g, seed_sur, 
  germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
  seed_disp_mat_1D, pollen_disp_mat, herb_app_reciv_expos);

herb_ee_ls_empty = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_empty, 
  source_Rr_reciv_empty, source_RR_reciv_empty, int_mean_g, int_sd_g, seed_sur, 
  germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
  seed_disp_mat_1D, pollen_disp_mat, herb_app_all_expos);
  
# herb_nn_ls_full = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_full, 
#   source_Rr_reciv_full, source_RR_reciv_full, int_mean_g, int_sd_g, seed_sur, 
#   germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
#   seed_disp_mat_1D, pollen_disp_mat, herb_app_none_expos);
# 
herb_en_ls_full = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_full, 
  source_Rr_reciv_full, source_RR_reciv_full, int_mean_g, int_sd_g, seed_sur, 
  germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
  seed_disp_mat_1D, pollen_disp_mat, herb_app_source_expos);

herb_ne_ls_full = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_full, 
  source_Rr_reciv_full, source_RR_reciv_full, int_mean_g, int_sd_g, seed_sur, 
  germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
  seed_disp_mat_1D, pollen_disp_mat, herb_app_reciv_expos);

herb_ee_ls_full = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_full, 
  source_Rr_reciv_full, source_RR_reciv_full, int_mean_g, int_sd_g, seed_sur, 
  germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
  seed_disp_mat_1D, pollen_disp_mat, herb_app_all_expos);

  
###########################################################################################
## PLOT THE RESULTS
# start getting the totnum for each scenario as will need this for all plots
# totnum_nn_empty = totnum_time_space(herb_nn_ls_empty, dg);
# totnum_nn_full = totnum_time_space(herb_nn_ls_full, dg);
totnum_ne_empty = totnum_time_space(herb_ne_ls_empty, dg);
totnum_ne_full = totnum_time_space(herb_ne_ls_full, dg);
totnum_en_empty = totnum_time_space(herb_en_ls_empty, dg);
totnum_en_full = totnum_time_space(herb_en_ls_full, dg);
totnum_ee_empty = totnum_time_space(herb_ee_ls_empty, dg);
totnum_ee_full = totnum_time_space(herb_ee_ls_full, dg);

totnum_max = maximum(hcat(totnum_ne_empty, totnum_ne_full, totnum_en_empty, totnum_en_full, 
  totnum_ee_empty, totnum_ee_full));
  
# pro_R summary
# proR_nn_empty = proR_time_space(herb_nn_ls_empty, dg);
# proR_nn_full = proR_time_space(herb_nn_ls_full, dg);
proR_ne_empty = proR_time_space(herb_ne_ls_empty, dg);
proR_ne_full = proR_time_space(herb_ne_ls_full, dg);
proR_en_empty = proR_time_space(herb_en_ls_empty, dg);
proR_en_full = proR_time_space(herb_en_ls_full, dg);
proR_ee_empty = proR_time_space(herb_ee_ls_empty, dg);
proR_ee_full = proR_time_space(herb_ee_ls_full, dg);

proR_max = maximum(hcat(proR_ne_empty, proR_ne_full, proR_en_empty, proR_en_full, proR_ee_empty, proR_ee_full));

# sur_rr_summary  
# surg_nn_empty = sur_g_time_space(herb_nn_ls_empty[3], dg, g_vals, herb_effect, 
#   base_sur, g_prot);
# surg_nn_full = sur_g_time_space(herb_nn_ls_full[3], dg, g_vals, herb_effect, 
#   base_sur, g_prot);
surg_ne_empty = sur_g_time_space(herb_ne_ls_empty[3], dg, g_vals, herb_effect, 
  base_sur, g_prot);
surg_ne_full = sur_g_time_space(herb_ne_ls_full[3], dg, g_vals, herb_effect, 
  base_sur, g_prot);
surg_en_empty = sur_g_time_space(herb_en_ls_empty[3], dg, g_vals, herb_effect, 
  base_sur, g_prot);
surg_en_full = sur_g_time_space(herb_en_ls_full[3], dg, g_vals, herb_effect, 
  base_sur, g_prot);
surg_ee_empty = sur_g_time_space(herb_ee_ls_empty[3], dg, g_vals, herb_effect, 
  base_sur, g_prot);
surg_ee_full = sur_g_time_space(herb_ee_ls_full[3], dg, g_vals, herb_effect, 
  base_sur, g_prot);

# TSR realative advantage, aslo create a mask of 0's and then switch the 0's to 1's
# relying on the mask to make the empty regions white
# TSR_adv_nn_empty = TSR_adv_time_space(herb_nn_ls_empty, dg, g_vals, germ_prob, fec_max,
#   g_effect_fec, dd_fec, sur_tup);
# mask_nn_empty = ones(size(TSR_adv_nn_empty));
# mask_nn_empty[TSR_adv_nn_empty .== 0.0] = 0.0;
# TSR_adv_nn_empty[TSR_adv_nn_empty .== 0.0] = 1.0;
# 
# TSR_adv_nn_full = TSR_adv_time_space(herb_nn_ls_full, dg, g_vals, germ_prob, fec_max,
#   g_effect_fec, dd_fec, sur_tup);
# mask_nn_full = ones(size(TSR_adv_nn_full));
# mask_nn_full[TSR_adv_nn_full .== 0.0] = 0.0;
# TSR_adv_nn_full[TSR_adv_nn_full .== 0.0] = 1.0;
# 
TSR_adv_ne_empty = TSR_adv_time_space(herb_ne_ls_empty, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec, sur_tup);
mask_ne_empty = ones(size(TSR_adv_ne_empty));
mask_ne_empty[TSR_adv_ne_empty .== 0.0] = 0.0;
TSR_adv_ne_empty[TSR_adv_ne_empty .== 0.0] = 1.0;

TSR_adv_ne_full = TSR_adv_time_space(herb_ne_ls_full, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec, sur_tup);
mask_ne_full = ones(size(TSR_adv_ne_full));
mask_ne_full[TSR_adv_ne_full .== 0.0] = 0.0;
TSR_adv_ne_full[TSR_adv_ne_full .== 0.0] = 1.0;

TSR_adv_en_empty = TSR_adv_time_space(herb_en_ls_empty, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec, sur_tup);
mask_en_empty = ones(size(TSR_adv_en_empty));
mask_en_empty[TSR_adv_en_empty .== 0.0] = 0.0;
TSR_adv_en_empty[TSR_adv_en_empty .== 0.0] = 1.0;

TSR_adv_en_full = TSR_adv_time_space(herb_en_ls_full, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec, sur_tup);
mask_en_full = ones(size(TSR_adv_en_full));
mask_en_full[TSR_adv_en_full .== 0.0] = 0.0;
TSR_adv_en_full[TSR_adv_en_full .== 0.0] = 1.0;

TSR_adv_ee_empty = TSR_adv_time_space(herb_ee_ls_empty, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec, sur_tup);
mask_ee_empty = ones(size(TSR_adv_ee_empty));
mask_ee_empty[TSR_adv_ee_empty .== 0.0] = 0.0;
TSR_adv_ee_empty[TSR_adv_ee_empty .== 0.0] = 1.0;

TSR_adv_ee_full = TSR_adv_time_space(herb_ee_ls_full, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec, sur_tup);
mask_ee_full = ones(size(TSR_adv_ee_full));
mask_ee_full[TSR_adv_ee_full .== 0.0] = 0.0;
TSR_adv_ee_full[TSR_adv_ee_full .== 0.0] = 1.0;
  
TSR_adv_max = maximum(hcat(TSR_adv_ne_empty, TSR_adv_ne_full, 
  TSR_adv_en_full, TSR_adv_ee_empty, TSR_adv_ee_full));
  
TSR_adv_min = minimum(hcat(TSR_adv_ne_empty, TSR_adv_ne_full, 
  TSR_adv_en_full, TSR_adv_ee_empty, TSR_adv_ee_full));
  
# TSR fec advantage, aslo create a mask of 0's and then switch the 0's to 1's
# relying on the mask to make the empty regions white
# TSR_rel_fec_nn_empty = TSR_rel_fec_time_space(herb_nn_ls_empty, dg, g_vals, germ_prob, fec_max,
#   g_effect_fec, dd_fec);
# maskrf_nn_empty = ones(size(TSR_rel_fec_nn_empty));
# maskrf_nn_empty[TSR_rel_fec_nn_empty .== 0.0] = 0.0;
# TSR_rel_fec_nn_empty[TSR_rel_fec_nn_empty .== 0.0] = 1.0;
# 
# TSR_rel_fec_nn_full = TSR_rel_fec_time_space(herb_nn_ls_full, dg, g_vals, germ_prob, fec_max,
#   g_effect_fec, dd_fec);
# maskrf_nn_full = ones(size(TSR_rel_fec_nn_full));
# maskrf_nn_full[TSR_rel_fec_nn_full .== 0.0] = 0.0;
# TSR_rel_fec_nn_full[TSR_rel_fec_nn_full .== 0.0] = 1.0;
# 
TSR_rel_fec_ne_empty = TSR_rel_fec_time_space(herb_ne_ls_empty, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec);
maskrf_ne_empty = ones(size(TSR_rel_fec_ne_empty));
maskrf_ne_empty[TSR_rel_fec_ne_empty .== 0.0] = 0.0;
TSR_rel_fec_ne_empty[TSR_rel_fec_ne_empty .== 0.0] = 1.0;

TSR_rel_fec_ne_full = TSR_rel_fec_time_space(herb_ne_ls_full, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec);
maskrf_ne_full = ones(size(TSR_rel_fec_ne_full));
maskrf_ne_full[TSR_rel_fec_ne_full .== 0.0] = 0.0;
TSR_rel_fec_ne_full[TSR_rel_fec_ne_full .== 0.0] = 1.0;

TSR_rel_fec_en_empty = TSR_rel_fec_time_space(herb_en_ls_empty, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec);
maskrf_en_empty = ones(size(TSR_rel_fec_en_empty));
maskrf_en_empty[TSR_rel_fec_en_empty .== 0.0] = 0.0;
TSR_rel_fec_en_empty[TSR_rel_fec_en_empty .== 0.0] = 1.0;

TSR_rel_fec_en_full = TSR_rel_fec_time_space(herb_en_ls_full, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec);
maskrf_en_full = ones(size(TSR_rel_fec_en_full));
maskrf_en_full[TSR_rel_fec_en_full .== 0.0] = 0.0;
TSR_rel_fec_en_full[TSR_rel_fec_en_full .== 0.0] = 1.0;

TSR_rel_fec_ee_empty = TSR_rel_fec_time_space(herb_ee_ls_empty, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec);
maskrf_ee_empty = ones(size(TSR_rel_fec_ee_empty));
maskrf_ee_empty[TSR_rel_fec_ee_empty .== 0.0] = 0.0;
TSR_rel_fec_ee_empty[TSR_rel_fec_ee_empty .== 0.0] = 1.0;

TSR_rel_fec_ee_full = TSR_rel_fec_time_space(herb_ee_ls_full, dg, g_vals, germ_prob, fec_max,
  g_effect_fec, dd_fec);
maskrf_ee_full = ones(size(TSR_rel_fec_ee_full));
maskrf_ee_full[TSR_rel_fec_ee_full .== 0.0] = 0.0;
TSR_rel_fec_ee_full[TSR_rel_fec_ee_full .== 0.0] = 1.0;
  
TSR_rel_fec_max = maximum(hcat(TSR_rel_fec_ne_empty, TSR_rel_fec_ne_full,  
  TSR_rel_fec_en_full, TSR_rel_fec_ee_empty, TSR_rel_fec_ee_full));
  
TSR_rel_fec_min = minimum(hcat(TSR_rel_fec_ne_empty, TSR_rel_fec_ne_full,
  TSR_rel_fec_en_full, TSR_rel_fec_ee_empty, TSR_rel_fec_ee_full));
#
###################################################################################################
# make the plots
  
# start by making the colormats for pro_R
colmat_ne_empty = colmat_2channel(totnum_ne_empty, proR_ne_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ne_full = colmat_2channel(totnum_ne_full, proR_ne_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_en_empty = colmat_2channel(totnum_en_empty, proR_en_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_en_full = colmat_2channel(totnum_en_full, proR_en_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ee_empty = colmat_2channel(totnum_ee_empty, proR_ee_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ee_full = colmat_2channel(totnum_ee_full, proR_ee_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);

  
dualchan_heatmap_grid(colmat_ee_empty, colmat_ee_full, colmat_ne_empty, colmat_ne_full, 
  colmat_en_empty, colmat_en_full, 1.9, totnum_max, 0.0, 1.0, "%R", 175.05, 360.0, output_loc, 
  "pro_R_time_space.pdf")

# survival colormats
colmat_ne_empty = colmat_2channel(totnum_ne_empty, surg_ne_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ne_full = colmat_2channel(totnum_ne_full, surg_ne_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_en_empty = colmat_2channel(totnum_en_empty, surg_en_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_en_full = colmat_2channel(totnum_en_full, surg_en_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ee_empty = colmat_2channel(totnum_ee_empty, surg_ee_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ee_full = colmat_2channel(totnum_ee_full, surg_ee_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);

dualchan_heatmap_grid(colmat_ee_empty, colmat_ee_full, colmat_ne_empty, colmat_ne_full, 
  colmat_en_empty, colmat_en_full, 1.9, totnum_max, 0.0, 1.0, "survival rr", 175.05, 360.0, output_loc, 
  "sur_rr_time_space.pdf")
 
# survival_rr and %R interaction colormats
colmat_ne_empty = colmat_2channel(totnum_ne_empty, proR_ne_empty .* surg_ne_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ne_full = colmat_2channel(totnum_ne_full,  proR_ne_full .* surg_ne_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_en_empty = colmat_2channel(totnum_en_empty,  proR_en_empty .* surg_en_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_en_full = colmat_2channel(totnum_en_full,  proR_en_full .* surg_en_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ee_empty = colmat_2channel(totnum_ee_empty,  proR_ee_empty .* surg_ee_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_ee_full = colmat_2channel(totnum_ee_full,  proR_ee_full .* surg_ee_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);

dualchan_heatmap_grid(colmat_ee_empty, colmat_ee_full, colmat_ne_empty, colmat_ne_full, 
  colmat_en_empty, colmat_en_full, 1.9, totnum_max, 0.0, 1.0, "survival rr x %R", 175.05, 360.0, output_loc, 
  "sur_rr_and_proR_time_space.pdf")

# TSR reproductive advantage
colmat_ne_empty = colmat_2channel(mask_ne_empty, TSR_adv_ne_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_adv_min, TSR_adv_max);
colmat_ne_full = colmat_2channel(mask_ne_full, TSR_adv_ne_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_adv_min, TSR_adv_max);
colmat_en_empty = colmat_2channel(mask_en_empty, TSR_adv_en_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_adv_min, TSR_adv_max);
colmat_en_full = colmat_2channel(mask_en_full, TSR_adv_en_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_adv_min, TSR_adv_max);
colmat_ee_empty = colmat_2channel(mask_ee_empty, TSR_adv_ee_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_adv_min, TSR_adv_max);
colmat_ee_full = colmat_2channel(mask_ee_full, TSR_adv_ee_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_adv_min, TSR_adv_max);
  
dualchan_heatmap_grid(colmat_ee_empty, colmat_ee_full, colmat_ne_empty, colmat_ne_full, 
  colmat_en_empty, colmat_en_full, 1.9, 1.0, TSR_adv_min, TSR_adv_max, "TSR reproductive advantage (K)", 
  175.05, 360.0, output_loc, "TSR_adv_time_space.pdf")
  
  
# TSR realative fecundity
colmat_ne_empty = colmat_2channel(maskrf_ne_empty, TSR_rel_fec_ne_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_rel_fec_min, TSR_rel_fec_max);
colmat_ne_full = colmat_2channel(maskrf_ne_full, TSR_rel_fec_ne_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_rel_fec_min, TSR_rel_fec_max);
colmat_en_empty = colmat_2channel(maskrf_en_empty, TSR_rel_fec_en_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_rel_fec_min, TSR_rel_fec_max);
colmat_en_full = colmat_2channel(maskrf_en_full, TSR_rel_fec_en_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_rel_fec_min, TSR_rel_fec_max);
colmat_ee_empty = colmat_2channel(maskrf_ee_empty, TSR_rel_fec_ee_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_rel_fec_min, TSR_rel_fec_max);
colmat_ee_full = colmat_2channel(maskrf_ee_full, TSR_rel_fec_ee_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, 1.0, TSR_rel_fec_min, TSR_rel_fec_max);
  
dualchan_heatmap_grid(colmat_ee_empty, colmat_ee_full, colmat_ne_empty, colmat_ne_full, 
  colmat_en_empty, colmat_en_full, 1.9, 1.0, TSR_rel_fec_min, TSR_rel_fec_max, "TSR fec / TSS fec", 
  175.05, 360.0, output_loc, "TSR_relfec_time_space.pdf")
  
  
# want to show the effect of both fec_cost and g_prot as they both interact
df_gpro_fr = DataFrame(g_pros = repeat(collect(0 : 0.125 : 2.5), inner = 21),
  fr = repeat(collect(0 : 0.05 : 1.0), outer = 21), proR_AUS_full = 0.0, 
  TSR_quant_AUS_full = 0.0, proR_AUS_empty = 0.0, TSR_quant_AUS_empty = 0.0);

# do the parameter tests and record the AUS indexes  
for i in 1:size(df_gpro_fr)[1]
 
  g_effect_fec = exp(-(fec0 - abs(g_vals) * df_gpro_fr[:fr][i]));
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, df_gpro_fr[:g_pros][i], pro_exposed);
  
  hold_ee_full = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_full, 
    source_Rr_reciv_full, source_RR_reciv_full, int_mean_g, int_sd_g, seed_sur, 
    germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
    seed_disp_mat_1D, pollen_disp_mat, herb_app_all_expos);

  hold_ee_empty = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_empty, 
    source_Rr_reciv_empty, source_RR_reciv_empty, int_mean_g, int_sd_g, seed_sur, 
    germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
    seed_disp_mat_1D, pollen_disp_mat, herb_app_all_expos);
    
  df_gpro_fr[:TSR_quant_AUS_full][i] = AUS_proR_surg(hold_ee_full, dg, dx,
    g_vals, herb_effect, base_sur, df_gpro_fr[:g_pros][i]);     
    
  df_gpro_fr[:proR_AUS_full][i] = AUS_proR(hold_ee_full, dg, dx);
 
  df_gpro_fr[:TSR_quant_AUS_empty][i] = AUS_proR_surg(hold_ee_empty, dg, dx,
    g_vals, herb_effect, base_sur, df_gpro_fr[:g_pros][i]);     
    
  df_gpro_fr[:proR_AUS_empty][i] = AUS_proR(hold_ee_empty, dg, dx);
  
  print(i)
  print("\n")
  
end


# This dataframe takes quiet a bit to produce, save results so I don't need to keep re-runing the model_output
cd(output_loc)
#writetable("gpro_fr_sweep.csv", df_gpro_fr, header = true)
df_gpro_fr = readtable("gpro_fr_sweep.csv", header = true)

par_sweep_heatmap(df_gpro_fr, 1.3, :g_pros, convert(String, L"$\rho$"), :fr, convert(String, L"$f_r$"), [:proR_AUS_empty, :proR_AUS_full, 
  :TSR_quant_AUS_empty, :TSR_quant_AUS_full], 0.0, 1.0, ["AUC %R", "AUC %R x sur rr"], 175.05, 360.0, 
  output_loc, "g_pro_fr_par_sweep.pdf")

  
# grid of two channel heat maps to show how populations change over time
# pass in 6 colormatircies to plot 

colmat_pR_empty = colmat_2channel(totnum_ee_empty, proR_ee_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_pR_full = colmat_2channel(totnum_ee_full, proR_ee_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_srr_empty = colmat_2channel(totnum_ee_empty, surg_ee_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_srr_full = colmat_2channel(totnum_ee_full, surg_ee_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_both_empty = colmat_2channel(totnum_ee_empty,  proR_ee_empty .* surg_ee_empty, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);
colmat_both_full = colmat_2channel(totnum_ee_full,  proR_ee_full .* surg_ee_full, 1.0, 0.5, 175.05, 360.0, 
  0.0, totnum_max, 0.0, 1.0);

dualchan_heatmap_3measure(colmat_pR_empty, colmat_pR_full,
  colmat_srr_empty, colmat_srr_full, colmat_both_empty,
  colmat_both_full, 1.9, 1.0, 0.0, 1.0, 
  175.05, 360.0, output_loc, "exposed_empty_3metrics.pdf")
  