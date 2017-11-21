# main calling script for TSR optimisation 

using HDF5, JLD

##########################################################################

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"
data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out"

cd(code_loc)
include("pop_process_2sb_2herb.jl"); 
include("managment_functions.jl"); 

##########################################################################
# inital populations 
int_SB1 = 100.0;
int_SB2 = 100.0;

int_RR = 0.0;
int_Rr = 0.01;
int_AA = 0.0;
int_Aa = 0.01;

# define parameters 
pro_exposed = 0.8;
sur_base = 10.0;
sur_h1 = 0.01; 
sur_h2 = 0.01;
sur_crop_alt = 0.8;
sur_spot = 0.05;
T = 20;
inv_frac = 0.8;
germ_prob = 0.5;
seed_sur = 0.5;
fec_max = 60.0; 
fec_dd = 0.001;

# reward function parameters
dis_rate = 0.96;
Y0 =1668.0;
Y_slope = 0.0032593;
Y_ALT = 769.0; # based on spring barley NIX 2017 
rep_pen = 0.92; # based on NIX 2017 pp. 9 point 7 

#cost parameters
cost_herb_one = 96.0; # cost herbicide from NIX 2017
cost_FAL = 36.0; # based on two glyphosate applications 
cost_WW = 383.0;
cost_ALT = 273.0 # based on costs from spring barley NIX 2017
cost_spot = 0.03; # clearing 100000 would cost £2000, > yeild of WW
cost_plow = 74.0; # cost inversion plowing

# set up the populations and mixing kernel
N_G = 9;

mk = make_TSR_mix_index(N_G);

pop = make_int_pop(int_SB1, int_SB2, int_RR, int_Rr, int_AA, int_Aa, mk);
offspring = zeros(length(pop[1]) ^ 2);

mix_kern = make_TSR_kernel(mk);


@time seeds = TSR_cross_ind(pop[1], pop[1] / sum(pop[1]), offspring, mk);

@time TSR_par_cross(pop[1], pop[1] / sum(pop[1]), offspring, mk);

@time mix_seed = mix_kern * offspring

