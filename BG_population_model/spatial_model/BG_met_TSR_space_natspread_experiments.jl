# natural spread experiments 
#######################################################################################################
# make a wrapper to run the natural spread scenarios in parallel  

########################################################################################################3
using DataFrames

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

x_dim = 500; # number of spatial evaluation points, actual landscape size is x_dim * dx
dx = 1.0;
source_locs = 1:20 # source population locations, at edge of landscape
recive_locs = 21:x_dim

int_num_RR = 0.1;
int_num_Rr = 0.0;
int_num_rr = 10.0; # number of intial seeds at each location for each genoptype, assume only TS susceptible
num_iter = 50;

seed_pro_short = 0.4; 
seed_mean_dist_short = 0.5; 
pro_seeds_to_mean_short = 0.4; 
seed_mean_dist_long = 1.5; 
pro_seeds_to_mean_long = 0.4;
scale_pollen = 32.0;
shape_pollen = 3.23; 
offspring_sd = 1.0;
fec_max = 45.0;
fec0 = 10.0;
fec_cost = 0.3;
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



run_natspread(g_vals::Array{Float64, 1}, x_dim::Int64, dg::Float64, 
  num_iter::Int64, int_rr::Array{Float64, 1}, int_Rr::Array{Float64, 1}, 
  int_RR::Array{Float64, 1}, int_mean_g::Float64, int_sd_g::Float64, seed_sur::Float64, 
  germ_prob::Float64, resist_G::Array{String, 1}, fec_max::Float64, dd_fec::Float64, 
  g_effect_fec::Array{Float64, 1}, sur_tup::Tuple{Float64, Array{Float64, 1}}, 
  offspring_sd::Float64, seed_disp_mat_1D::Array{Float64, 2}, pollen_disp_mat::Array{Float64, 2}, 
  herb_app::Array{Int64, 1})
#Little test example down here to figure out how to plot a 2 color channle matrix 

# heat map attempt
heatmap([RGBA(0.0, 0.0, 0.0, 1.0) RGBA(0.4, 1.0, 0.0, 1.0); 
  RGBA(1.0, 0.4, 0.0, 1.0) RGBA(0.0, 0.0, 3.0, 1.0)])


