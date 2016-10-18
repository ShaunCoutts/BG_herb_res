# simulation experiments for the TSR and NTSR joint evolution
# to test how TSR can invade into a population with NTSR

# Ideally I will be be able to specify which parameters to vary and
# the values to test them at. I could just pass every parameter 
# then I get length and run a for loop over it.

# need to run each source scenario three times, once in a empty landscape, one in an exposed population
# and one in a herbicide exposed population (can pre calculate all these so don't need to re-caclulate each time 
file_loc_func_p = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model" 
cd(file_loc_func_p)
include("BG_met_TSR_space_exper_funs.jl")
include("BG_met_TSR_space_pop_process.jl")
include("BG_met_TSR_space_dispersal_functions.jl")
  


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

int_num_RR = 0.0;
int_num_Rr = 0.0;
int_num_rr = 10.0; # number of intial seeds at each location for each genoptype, assume only TS susceptible
burnin = 20;
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
fec_cost = 1.0;
dd_fec = 0.025;
base_sur = 10.0; 
herb_effect = 20.0; 
g_prot = 1.0; 
pro_exposed = 0.8;
seed_sur = 0.45;
germ_prob = 0.52;
resist_G = ["RR", "Rr"];

# set up the evaluation points for quantitative resistance
g_vals = collect(lower_g : dg : upper_g);   


# The set up the set of seeds added after the burnin period
# Four different source populations: low mean_g and low %R, low mean_g and high %R, 
# high mean_g and low %R, high mean_g and high %R

# low mean_g = 0 
inject_g_low = 0.0;
# high mean_g = g such that s = 0.95 for exposed individuals 
inject_g_high = sur_2_g(0.95, herb_effect, base_sur, g_prot); 

inject_sd_g = int_sd_g;

#low and high TSR
int_TSR_low = 0.05;
int_TSR_high = 0.95;

num_inject = 10.0;
inject_locs = [convert(Int64, x_dim / 2)];

# Do the low g and low TSR scenario
g_low_TSR_low = run_scene_trans(g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, int_TSR_low,
  inject_g_low, inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
  resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, g_prot, pro_exposed, 
  seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
  pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd);










1 / (1 + exp(-(base_sur - (herb_effect - int_g_high * g_prot))))

run_scene_trans(g_vals::Array{Float64, 1}, x_dim::Int64, dg::Float64, dx::Float64, 
  num_iter::Int64, burnin::Int64, num_inject::Float64, pro_R_inject::Float64, inject_mean_g::Float64, 
  inject_sd_g::Float64, inject_locs::Array{Float64, 1}, int_rr::Float64, int_mean_g::Float64,
  int_sd_g::Float64, seed_sur::Float64, germ_prob::Float64, resist_G::Array{ASCIIString, 1}, 
  fec_max::Float64, dd_fec::Float64, fec0::Float64, fec_cost::Float64, base_sur::Float64, 
  herb_effect::Float64, g_prot::Float64, pro_exposed::Float64, seed_pro_short::Float64, 
  seed_mean_dist_shor::Float64, pro_seeds_to_mean_short::Float64, seed_mean_dist_long::Float64, 
  pro_seeds_to_mean_long::Float64, scale_pollen::Float64, shape_pollen::Float64, offspring_sd::Float64)














