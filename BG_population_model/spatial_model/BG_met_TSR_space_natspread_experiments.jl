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

seed_pro_short = 0.4; 
seed_mean_dist_short = 0.5; 
pro_seeds_to_mean_short = 0.4; 
seed_mean_dist_long = 5.0; 
pro_seeds_to_mean_long = 0.2;
scale_pollen = 32.0;
shape_pollen = 3.23; 
offspring_sd = 1.0;
fec_max = 60.0;
fec0 = 10.0;
fec_cost = 1.0;
dd_fec = 0.15;
base_sur = 10.0; 
herb_effect = 15.0; 
g_prot = 2.0; 
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
herb_nn_ls_empty = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_empty, 
  source_Rr_reciv_empty, source_RR_reciv_empty, int_mean_g, int_sd_g, seed_sur, 
  germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
  seed_disp_mat_1D, pollen_disp_mat, herb_app_none_expos);

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
  
herb_nn_ls_full = run_natspread(g_vals, x_dim, dg, num_iter, source_rr_reciv_full, 
  source_Rr_reciv_full, source_RR_reciv_full, int_mean_g, int_sd_g, seed_sur, 
  germ_prob, resist_G, fec_max, dd_fec, g_effect_fec, sur_tup, offspring_sd, 
  seed_disp_mat_1D, pollen_disp_mat, herb_app_none_expos);

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
totnum_nn_empty = totnum_time_space(herb_nn_ls_empty, dg);
totnum_nn_full = totnum_time_space(herb_nn_ls_full, dg);
totnum_ne_empty = totnum_time_space(herb_ne_ls_empty, dg);
totnum_ne_full = totnum_time_space(herb_ne_ls_full, dg);
totnum_en_empty = totnum_time_space(herb_en_ls_empty, dg);
totnum_en_full = totnum_time_space(herb_en_ls_full, dg);
totnum_ee_empty = totnum_time_space(herb_ee_ls_empty, dg);
totnum_ee_full = totnum_time_space(herb_ee_ls_full, dg);

totnum_max = maximum(hcat(totnum_nn_empty, totnum_nn_full, totnum_ne_empty,
  totnum_ne_full, totnum_en_empty, totnum_en_full, totnum_ee_empty, totnum_ee_full));
  
# pro_R summary
proR_nn_empty = proR_time_space(herb_nn_ls_empty, dg);
proR_nn_full = proR_time_space(herb_nn_ls_full, dg);
proR_ne_empty = proR_time_space(herb_ne_ls_empty, dg);
proR_ne_full = proR_time_space(herb_ne_ls_full, dg);
proR_en_empty = proR_time_space(herb_en_ls_empty, dg);
proR_en_full = proR_time_space(herb_en_ls_full, dg);
proR_ee_empty = proR_time_space(herb_ee_ls_empty, dg);
proR_ee_full = proR_time_space(herb_ee_ls_full, dg);

proR_max = maximum(hcat(proR_nn_empty, proR_nn_full, proR_ne_empty,
  proR_ne_full, proR_en_empty, proR_en_full, proR_ee_empty, proR_ee_full));

# sur_rr_summary  
surg_nn_empty = sur_g_time_space(herb_nn_ls_empty[3], dg, g_vals, herb_effect, 
  base_sur, g_prot);
surg_nn_full = sur_g_time_space(herb_nn_ls_full[3], dg, g_vals, herb_effect, 
  base_sur, g_prot);
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
 