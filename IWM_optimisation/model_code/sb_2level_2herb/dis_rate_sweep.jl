# parameter sweep done in parallel to spped things up as solving for each 
# parameter combination will probably take around 4hrs, so we will have 
# be strategic around what we test

# main calling script 
using HDF5, JLD
#@everywhere using StatsBase

# file locations 
#@everywhere code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"
#@everywhere data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out"

# file locations big mac
@everywhere code_loc = "/Users/shauncoutts/IWM_optim/model_code/"
@everywhere out_loc = "/Users/shauncoutts/IWM_optim/out_obj/"

# file locations ICEBERG
#@everywhere code_loc = "/home/bo1src/IWM_optim/model_code/"
#@everywhere out_loc = "/home/bo1src/IWM_optim/output/"
#
##########################################################################
# get reuired helper functions
@everywhere cd(code_loc)
@everywhere include("pop_process_2sb_2herb.jl"); 
@everywhere include("managment_functions.jl"); 

##########################################################################
# MAKE A WRAPPER FUNCTION THAT TAKES A LIST OF PARAMETERS AND 
# RUNS THE SOLVER, RETURNING THE TUPLE SOLITION. 
# Send parameters as dataframe so the names are easier to keep track of
@everywhere function GA_wrapper(p_tup::Tuple{Dict{Symbol, Any}, Int64})

	p = p_tup[1]
	# create a string to print out to see progress
	iter = p_tup[2]
	print("$iter \n")

	sol = GA_solve(p[:T], p[:pop_size], p[:num_gen], 
		p[:cost_herb_one], p[:cost_WW], p[:cost_ALT], 
		p[:cost_FAL], p[:cost_plow], p[:spot_fix], p[:spot_var], 
		p[:sur_crop_alt], p[:low_g], p[:up_g], p[:dg], 
		p[:off_sd], p[:off_cv], p[:int_N], p[:int_sd], 
		p[:int_cv], p[:int_g1], p[:int_g2], p[:inv_frac], 
		p[:germ_prob], p[:seed_sur], p[:fec_max], p[:fec_dd], 
	       	p[:sur_spot], p[:dis_rate], p[:Y0], p[:Y_slope],
		p[:Y_ALT], p[:pro_exposed], p[:sur_base], p[:rep_pen], 
		p[:effect_herb1], p[:effect_herb2], p[:prot_g1_herb1], 
		p[:prot_g2_herb2], p[:fr], p[:f0], p[:mut]);

	return sol

end

##########################################################################
# fixed parameters
# define resistance trait values and their co-var
low_g = -10.0;
dg = 1.0;
up_g = 10.0;
off_sd = 1.1225; # Va = 1.5

# define parameters 
pro_exposed = 0.8;
sur_base = 10.0;
effect_herb1 = 14.0; 
effect_herb2 = 14.0;
prot_g1_herb1 = 1.3; # these have to be in a pretty narrow range to get 
prot_g2_herb2 = 1.3; # sensible results. It all interacts with fr and f0.
sur_crop_alt = 0.8;
sur_spot = 0.05;
T = 20;
inv_frac = 0.8;
germ_prob = 0.5;
seed_sur = 0.5;
fec_max = 60.0; 
fec_dd = 0.00001;
fr = 0.2;
f0 = 4.0;

# reward function parameters
Y_slope = 0.0032593;
Y_ALT = 769.0; # based on spring barley NIX 2017 
rep_pen = 0.92; # based on NIX 2017 pp. 9 point 7 

#cost parameters
cost_herb_one = 96.0; # cost herbicide from NIX 2017
cost_FAL = 36.0; # based on two glyphosate applications 
cost_WW = 383.0;
cost_ALT = 273.0; # based on costs from spring barley NIX 2017
spot_fix = 20.0; # based on 1 ha weed contractor rates for non-boom
spot_var = 0.03; # clearing 100000 would cost £3000, > yeild of WW
cost_plow = 74.0; # cost inversion plowing

# initial conditions 
int_sd = 1.732;
int_cv = 0.0;

##########################################################################
# varied parameters as lists
off_cv = [0.0, 0.5];
dis_rate = [0.85, 0.9, 0.95, 1.0];
Y0 = [1022, 1668.0];
int_N = [100000.0];
int_g1 = [sur_2_g(0.02, effect_herb1, sur_base, prot_g1_herb1),
	  sur_2_g(0.97, effect_herb1, sur_base, prot_g1_herb1),
	  sur_2_g(0.97, effect_herb1, sur_base, prot_g1_herb1)];
int_g2 = [sur_2_g(0.02, effect_herb2, sur_base, prot_g2_herb2),
	  sur_2_g(0.02, effect_herb2, sur_base, prot_g2_herb2),
	  sur_2_g(0.97, effect_herb2, sur_base, prot_g2_herb2)];

# GA parameters
num_gen = 70;
pop_size = 200; # num of action sequences evaluaed each gen, must be even
mut = 0.02;

# make an array of dataframes to pass to the wrapper function
par_dfs = [];
count = 1;
for dr in dis_rate
	for ocv in off_cv
		for y0 in Y0
			for int_n in int_N 
				for ig in 1:length(int_g1)

				push!(par_dfs, 
				      	(Dict(:T => T, 
					      :pop_size => pop_size, 
					      :num_gen => num_gen, 
					      :mut => mut, 
					      :cost_WW => cost_WW, 
					      :cost_ALT => cost_ALT, 
					      :cost_FAL => cost_FAL, 
					      :cost_herb_one => cost_herb_one, 
					      :cost_plow => cost_plow, 
					      :spot_fix => spot_fix,
					      :spot_var => spot_var,
					      :sur_crop_alt => sur_crop_alt, 
					      :low_g => low_g, 
					      :up_g => up_g, 
					      :dg => dg, 
					      :off_sd => off_sd, 
					      :off_cv => ocv, 
					      :int_N => int_n, 
					      :int_sd => int_sd, 
					      :int_cv => int_cv,
					      :int_g1 => int_g1[ig],
					      :int_g2 => int_g2[ig],
					      :inv_frac => inv_frac,
					      :germ_prob => germ_prob,
					      :seed_sur => seed_sur, 
					      :fec_max => fec_max, 
					      :fec_dd => fec_dd, 
					      :sur_spot => sur_spot, 
					      :dis_rate => dr, 
					      :Y0 => y0, 
					      :Y_slope => Y_slope, 
					      :Y_ALT => Y_ALT,
					      :pro_exposed => pro_exposed, 
					      :sur_base => sur_base,
					      :rep_pen => rep_pen, 
					      :effect_herb1 => effect_herb1,
					      :effect_herb2 => effect_herb2, 
					      :prot_g1_herb1 => prot_g1_herb1, 
					      :prot_g2_herb2 => prot_g2_herb2,
					      :fr => fr, 
					      :f0 => f0), 
					 count))

					count += 1;

				end
			end
		end
	end
end

dis_sweep = pmap(GA_wrapper, par_dfs)	

cd(out_loc)

save("dis_sweep_obj.jld", "dis_sweep", dis_sweep)

#sol_sweep = load("sol_sweep.jld")["sol_sweep"]

