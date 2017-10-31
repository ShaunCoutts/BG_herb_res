# main calling script 

##########################################################################

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"

cd(code_loc)
include("pop_process_2sb_2herb.jl"); 
include("managment_functions.jl"); 

# define resistance trait values and their co-var
low_g = -15.0;
dg = 1.5;
up_g = 15.0;
off_sd = 1.141;
off_cv = 0.0;

# define parameters 
pro_exposed = 0.8;
sur_base = 10.0;
effect_herb1 = 15.0; 
effect_herb2 = 15.0;
prot_g1_herb1 = 2.0; 
prot_g1_herb2 = 1.0;
prot_g2_herb1 = 1.0;
prot_g2_herb2 = 2.0;
sur_crop_alt = 0.8;
T = 20;
inv_frac = 0.8;
germ_prob = 0.5;
seed_sur = 0.5;
fec_max = 60.0; 
fec_dd = 0.001;
fr = 0.5;
f0 = 4.0;
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
cost_spot = 0.02; # clearing 100000 would cost £2000, > yeild of WW
cost_plow = 74.0; # cost inversion plowing

# initial conditions 
int_g1 = 0.0;
int_g2 = 0.0;
int_sd = 1.141;
int_cv = 0.0;
int_N = 10.0;

# GA parameters
num_gen = 30;
pop_size = 50; # population of action sequences 
mut = 0.02;

@time sol = GA_solve(T, pop_size, num_gen, cost_herb_one, cost_WW, cost_ALT, cost_FAL, 
	cost_plow, cost_spot, sur_crop_alt, low_g, up_g, dg, off_sd, off_cv, 
	int_N, int_sd, int_cv, int_g1, int_g2, inv_frac, germ_prob, seed_sur, 
	fec_max, fec_dd, dis_rate, Y0,Y_slope, Y_ALT, pro_exposed, sur_base, 
	rep_pen, effect_herb1, effect_herb2, prot_g1_herb1, prot_g1_herb2, 
	prot_g2_herb1, prot_g2_herb2, fr, f0, mut)











# define precomputed objects
crop_sur_tup = (1.0, sur_crop_alt, 0.0);
spot_sur_tup = (0.0, 1.0);

herb_sur_tup = survial_herb_setup(g1_vals, g2_vals, pro_exposed, sur_base, 
  effect_herb1, effect_herb2, prot_g1_herb1, prot_g1_herb2, prot_g2_herb1, prot_g2_herb2);

mix_keys = make_index_keys(len_g, len_g);

off_tem = offspring_dist_setup(g1_vals, g2_vals, cov_mat, mix_keys);

fec_cost = fec_cost_maker(fr, f0, g1_vals, g2_vals);

cost_space = make_cost_space(cost_herb_fixed, cost_herb_var, cost_fal, cost_plow);


# define holding arrays
seedbank_1 = zeros(time_horizion + 1, size(g1_vals)[1]);
seedbank_2 = zeros(time_horizion + 1, size(g1_vals)[1]);
mat_pop = zeros(time_horizion + 1, size(g1_vals)[1]);
pat_pop = zeros(size(g1_vals)[1]);
new_seeds = zeros(size(g1_vals)[1]);
par_mix = zeros(size(mix_keys[1]));

action_space = make_action_space();
act_seq = act_seq_herb0(action_space, time_horizion)
#make the plow array
plow_tup = plow_subact[act_seq[:, ACT_PLOW]];
plow_seq = collect(plow_tup);

# create an intial population 
int_cov = [1.414 0.0; 0.0 1.414];
int_dist = MvNormal(int_cov);
int_sb1 = pdf(int_dist, transpose(hcat(g1_vals, g2_vals))) 
int_sb2 = zeros(size(g1_vals)[1]) 

# create set of action sequences
int_pop = rand_act_pop(length(action_space), time_horizion, 100);
# turn the first action seqence into sub-actions
acts = act_seq_2_sub_act(action_space, int_pop[1, :]);
plow_tup = plow_subact[acts[:, ACT_PLOW]];
plow_seq = collect(plow_tup);

# test one
@time one_run!(seedbank_1, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup, acts[:, ACT_CROP], plow_seq, acts[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg)

@time one_step!(seedbank_1[2, :], seedbank_2[2, :], mat_pop[2, :], pat_pop, new_seeds, 
  par_mix, mix_keys, off_tem, fec_cost, crop_sur_tup, herb_sur_tup, acts[2, ACT_CROP], 
  plow_seq[2], acts[2, ACT_HERB], inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg)




R = reward_total(mat_pop, dis_rate, Y0, Y_slope, Y_shape, N_half_Y, alt_Y, social_scale, 
  social_shape, econ_wieght, acts, cost_space, econ_constraint, dg)

  
a = costs(acts, cost_space)

tot_pop = vcat(sum(mat_pop, 2)...) * dg * dg;

ER = economic_reward(tot_pop[2:end], acts[:, ACT_CROP], Y0, Y_slope, Y_shape, 
    N_half_Y, alt_Y)

prof = 



