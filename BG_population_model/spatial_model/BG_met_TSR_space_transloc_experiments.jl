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
int_g_mean = 0.0;
int_g_sd = 1.4142;

x_dim = 500; # number of spatial evaluation points, actual landscape size is x_dim * dx
dx = 1.0;

int_num_RR = 0;
int_num_Rr = 0;
int_num_rr = 10; # number of intial seeds at each location for each genoptype, assume only TS susceptible
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
fec0 = 10.0;
fec_cost = 1.0;
base_sur = 10.0; 
herb_effect = 20.0; 
g_prot = 1.0; 
pro_exposed = 0.8;


# set up the evaluation points for quantitative resistance
g_vals = collect(lower_g : dg : upper_g);   
















