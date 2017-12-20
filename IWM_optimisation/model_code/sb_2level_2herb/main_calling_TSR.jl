# main calling script for TSR optimisation 

using HDF5, JLD

##########################################################################

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"
data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out"

cd(code_loc)
include("pop_process_TSR.jl"); 
include("managment_functions.jl"); 

##########################################################################
# inital populations 
int_SB1 = 100.0;
int_SB2 = 100.0;

int_RR = 0.0;
int_Rr = 0.1;
int_AA = 0.0;
int_Aa = 0.01;
int_N = 100000.0;

# define parameters 
pr_ex = 0.8;
s0 = 0.99;
sur_herb = 0.01; 
sur_crop_alt = 0.2;
sur_spot = 0.05;
T = 20;
inv_frac = 0.6;
germ_prob = 0.5;
seed_sur = 0.5;
fec_max = 60.0; 
fec_dd = 0.0003;

# reward function parameters
dis_rate = 0.96;
Y0 =1668.0;
Y_slope = 0.0032593;
Y_ALT = 769.0; # based on spring barley NIX 2017 
rep_pen = 0.92; # based on NIX 2017 pp. 9 point 7 

#cost parameters
cost_herb = 96.0; # cost herbicide from NIX 2017
cost_FAL = 36.0; # based on two glyphosate applications 
cost_WW = 383.0;
cost_ALT = 273.0 # based on costs from spring barley NIX 2017
spot_fix = 20.0; # clearing 100000 would cost £2000, > yeild of WW
spot_var = 0.03; # clearing 100000 would cost £2000, > yeild of WW
cost_plow = 74.0; # cost inversion plowing

# set up the populations and mixing kernel
N_G = 9;
T = 20;

mk = make_TSR_mix_index(N_G);

pop = make_int_pop(int_SB1, int_SB2, int_RR, int_Rr, int_AA, int_Aa, mk);
offspring = zeros(length(pop[1]) ^ 2);

mix_kern = make_TSR_kernel(mk);

herb_sur = make_herb_sur_dom(mk, s0, pro_exposed, sur_h1);

@time TSR_par_cross!(pop[1], pop[1] / sum(pop[1]), offspring, mk);

@time mix_seed = mix_kern * offspring

# make some populations and holding arrays
SB1 = zeros(T + 1, N_G);
SB2 = zeros(T + 1, N_G);

SB1[1, :] = deepcopy(pop[1]);
SB2[1, :] = deepcopy(pop[2]);

SB1[2, :] = deepcopy(pop[1]);
SB2[2, :] = deepcopy(pop[2]);

ab_pop = zeros(T + 1, N_G);
ab_pop_spot = zeros(T + 1, N_G);

eff_pop = zeros(N_G);
mat_pop = zeros(N_G);
pat_pop = zeros(N_G);

# make some actions for testing 
herb_act = HERB1;
crop_act = CROP_WW;
plow = false;
spot = 0;

herb_act_seq = repeat([HERB12], inner = T);
crop_act_seq = repeat([CROP_WW], inner = T);
plow_seq = repeat([false], inner = T);
spot_act_seq = repeat([0], inner = T);

@time onestep!(SB1, SB2, ab_pop, ab_pop_spot, eff_pop, mat_pop, pat_pop, 
	offspring, mix_kern, mk, herb_sur, inv_frac, germ_prob,
	seed_sur, sur_crop_alt, sur_spot, fec_dd, fec_max,
	herb_act, crop_act, plow, spot, 2);

@time one_run!(SB1, SB2, ab_pop, ab_pop_spot, eff_pop, mat_pop, pat_pop, 
	offspring, mix_kern, mk, herb_sur, inv_frac, germ_prob, seed_sur, 
	sur_crop_alt, sur_spot, fec_dd, fec_max, herb_act_seq, crop_act_seq,
	plow_seq, spot_act_seq, pop[1], pop[2], T)
  

# test out the GA solver
pop_size = 1000;
num_gen = 100;
mut = 0.03;

@time sol =  GA_solve_TSR(T, pop_size, num_gen, mut,
	cost_herb, cost_WW, cost_ALT, cost_FAL, cost_plow, spot_fix, spot_var, 
	sur_crop_alt, sur_herb,
	int_N, int_RR, int_Rr, int_AA, int_Aa, 
	inv_frac, germ_prob, seed_sur, fec_max, 
	fec_dd, sur_spot, dis_rate, Y0,
	Y_slope,Y_ALT, pr_ex, s0, rep_pen);





