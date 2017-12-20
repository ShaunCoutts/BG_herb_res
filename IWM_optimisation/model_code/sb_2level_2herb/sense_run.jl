# main calling script for TSR parameter sweep where we cover alot of
# parameter space. Run times are managable on the 16 core machine but
# memory footprint of the output is very large if I hold the whole object
# would take 780 GB to do 50000 parameter combinations. Instead just record
# the parameters and the optimal action, then later can use this to simulate
# the population in the future to get perfomance metrics.

using HDF5, JLD
using BlackBoxOptim
using DataFrames

##########################################################################
# wrapper function to pass to the solver for parallization
@everywhere function GA_wrapper(p::Dict{Symbol, Any})

	sol = GA_solve_TSR(p[:T], p[:pop_size], p[:num_gen], p[:mut],
		p[:cost_herb], p[:cost_WW], p[:cost_ALT], p[:cost_FAL], 
		p[:cost_plow], p[:spot_fix], p[:spot_var], p[:sur_alt], 
		p[:sur_herb], p[:int_N], p[:int_RR], p[:int_Rr], 
		p[:int_AA], p[:int_Aa], p[:inv_frac], p[:germ_prob], 
		p[:seed_sur], p[:fec_max], p[:fec_dd], p[:sur_spot], 
		p[:dis_rate], p[:Y0], p[:Y_slope],p[:Y_ALT], p[:pr_ex], 
		p[:s0], p[:rep_pen])
	
	
	best_ind = best_n_seq(1, sol[2])
	out = merge(sol[3], Dict(:act_seq => sol[1][end][best_ind[end], :]))	

	return out

end

# helper function to logit transformation to help sample the proprotions
# in a more useful space 
function logit(p::Float64)

	return log(p) - log(1 - p)

end

function inv_logit(n::Float64)

	return 1 / (1 + exp(-n))

end

##########################################################################

#@everywhere code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"
#data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out"

@everywhere code_loc = "/Users/shauncoutts/IWM_optim/model_code/"
data_loc = "/Users/shauncoutts/IWM_optim/out_obj/"

@everywhere cd(code_loc)
@everywhere include("pop_process_TSR.jl"); 
@everywhere include("act_seq_handelers.jl"); 
@everywhere include("managment_functions.jl"); 

##########################################################################
# set the upper and lower limits of parameters to test over, hold as many
# constant as I can get away with, also think about the scales that are 
# being sampled on, to try and sample the most interesting space more intensly 
# inital populations 
int_RR = 0.0;
int_Rr_u_raw = 5.0;
int_Rr_l_raw = -10.0;

int_AA = 0.0;
int_Aa_u_raw = 5.0;
int_Aa_l_raw = -10.0;

int_N_u = 100000.0;
int_N_l = 100.0;

# define parameters 
pr_ex_u = 1.0; # unrealistically high
pr_ex_l = 0.7; # implies 30% of blackgrass is missed, seems high

s0 = 0.99; # keep this fixed and assume it is partially wraped 
	   # up in fecundity

sur_herb = 0.01; # how effective is herb on suceptable plants keep fixed

sur_alt_u = 1.0 - 0.96; # range from Lutman et al 2013
sur_alt_l = 1.0 - 0.78;

sur_spot_u = 0.2;
sur_spot_l = 0.05; # range not known so make  wide, know 0.05 is efective

inv_frac_u = 0.9; # upper from spader in Grundy et al 1999 transition mat
inv_frac_l = 0.5; # put the lower bound to be very ineffective 

germ_prob_u = 0.6; # values from observed values in Colbach et al 2006
germ_prob_l = 0.45; # values from observed values in Colbach et al 2006

seed_sur_u = 0.86; # the range for seed survival is one of the better 
seed_sur_l = 0.2;  # supported

fec_max_u = 300.0; # range from Moss power point 
fec_max_l = 30.0; 

fec_dd_u = 0.001;
fec_dd_l = 0.0001;

# reward function parameters
dis_rate_l = 0.75;
dis_rate_u = 1.0;

Y0_u = 1758.0; # upper is upper 95CI from data
Y0_l = 986.0; # lower is low output Nix 2017 value

Y_slope_u = 0.00620878; # range from 95CI fitted yield model
Y_slope_l = 0.0002117964; # range from 95CI fitted yield model

# Set this up based on gross margin, simply set the cost to 0.0
Y_ALT_u = 647.0; # based on spring barley NIX 2017 
Y_ALT_l = 399.0; # based on spring barley NIX 2017 

rep_pen_u = 1.0; # based on NIX 2017 pp. 9 point 7 
rep_pen_l = 0.85; 

#cost parameters
cost_herb_u = 100.0; # cost herbicide from NIX 2017
cost_herb_l = 50.0; # cost herbicide from NIX 2017

cost_FAL = 36.0; # based on two glyphosate applications 

cost_WW = 383.0; # from Nix 2017 excluding herbicide, no range let it be captued in Y0

cost_ALT = 0.0 # set to 0.0 so only sampling over gross margin

spot_fix_u = 100.0; # just set a wide range to see how it affect the decision 
spot_fix_l = 10.0;

spot_var_u = 0.05; # clearing 100000 would cost £5000, > yeild of WW
spot_var_l = 0.01; # clearing 100000 would cost £1000 

cost_plow_u = 74.0 * 1.25; # Nix 2017 pp 202 seems about a 23% chang in price 
cost_plow_l = 74.0 * 0.75; # between light and heavy soil, use estimate +- 0.25

# set up the populations and mixing kernel
N_G = 9;
T = 25;

pop_size = 1000;
num_gen = 100;
mut = 0.03;
num_par_combs = 50000;

# sample the parameter space with LHS
par0 = BlackBoxOptim.Utils.latin_hypercube_sampling(
	[int_Rr_l_raw, int_Aa_l_raw, int_N_l, pr_ex_l, sur_alt_l,
	 sur_spot_l, inv_frac_l, germ_prob_l, seed_sur_l, fec_max_l, 
	 fec_dd_l, dis_rate_l, Y0_l, Y_slope_l, Y_ALT_l, rep_pen_l,
	 cost_herb_l, spot_fix_l, spot_var_l, cost_plow_l], 
	[int_Rr_u_raw, int_Aa_u_raw, int_N_u, pr_ex_u, sur_alt_u,
	 sur_spot_u, inv_frac_u, germ_prob_u, seed_sur_u, fec_max_u, 
	 fec_dd_u, dis_rate_u, Y0_u, Y_slope_u, Y_ALT_u, rep_pen_u,
	 cost_herb_u, spot_fix_u, spot_var_u, cost_plow_u], 
	num_par_combs)

par_list = [];

for i in 1:num_par_combs

	push!(par_list,
	      Dict(:int_Rr => inv_logit(par0[1, i]),
		   :int_Aa => inv_logit(par0[2, i]),
		   :int_N => par0[3, i],
		   :pr_ex => par0[4, i],
		   :sur_alt => par0[5, i],
		   :sur_spot => par0[6, i],
		   :inv_frac => par0[7, i],
		   :germ_prob => par0[8, i],
		   :seed_sur => par0[9, i],
		   :fec_max => par0[10, i],
		   :fec_dd => par0[11, i],
		   :dis_rate => par0[12, i],
		   :Y0 => par0[13, i],
		   :Y_slope => par0[14, i],
		   :Y_ALT => par0[15, i],
		   :rep_pen => par0[16, i],
		   :cost_herb => par0[17, i],
		   :spot_fix => par0[18, i],
		   :spot_var => par0[19, i],
		   :cost_plow => par0[20, i],
		   :cost_WW => cost_WW,
		   :cost_ALT => cost_ALT,
		   :cost_FAL => cost_FAL,
		   :sur_herb => sur_herb,
		   :int_RR => int_RR,
		   :int_AA => int_AA,
		   :s0 => s0,
		   :T => T,
		   :pop_size => pop_size, 
		   :num_gen => num_gen,
		   :mut => mut))

end

@time sens_runs = pmap(GA_wrapper, par_list);

# save the resulting object to reuse later and it will be time consuming to create
cd(data_loc);
save("sens_runs_obj.jld", "sens_obj", sens_runs);

#sol_sweep = load("sens_runs_obj.jld")["sens_obj"]



