using DataFrames

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"

cd(code_loc)
include("pop_process_2sb_2herb.jl"); 
include("managment_functions.jl"); 

##########################################################################
# some helper function to parse and unpack the solver results

# takes an array of rewards, each element and array of rewards at each 
# generation. Return the indicies of n best seqences in each generation 
function best_n_seq(n::Int64, reward_arr::Array{Array, 1})

	N_gen = length(reward_arr)
	best_ind = Array{Int64, 2}(N_gen, n)

	for g in 1:N_gen

		si = sortperm(reward_arr[g], rev = true)
		best_ind[g, :] = si[1:n]
	
	end

	return best_ind

end

# pull out the rewards
function reward_gen(reward_arr::Array{Array, 1}, ind::Array{Int64, 2})

	N_gen = length(reward_arr)
	N_series = size(ind)[2]

	out = zeros(N_gen, N_series)
	# run through each series by generation, pulling out the rewards
	for g in 1:N_gen

		out[g, :] = reward_arr[g][ind[g, :]]
		
	end

	return out

end

# turn a seqence of action mumbers to seqences of herb, crop, plow, and spot
function actseq_subact(actseq::Array{Int64, 2}, A::Tuple)

	# holding arrays
	herb_seq = Array{Int64, 2}(size(actseq))
	crop_seq = Array{Int64, 2}(size(actseq))
	plow_seq = Array{Int64, 2}(size(actseq))
	spot_seq = Array{Int64, 2}(size(actseq))

	for s in 1:size(actseq)[1]
		
		acts = act_seq_2_sub_act(A, actseq[s, :])

		herb_seq[s, :] = acts[:, ACT_HERB]
		crop_seq[s, :] = acts[:, ACT_CROP]
		plow_seq[s, :] = acts[:, ACT_PLOW]
		spot_seq[s, :] = acts[:, ACT_SPOT]

	end

	return (herb_seq, crop_seq, plow_seq, spot_seq)

end

# take the whole solution and output the best performing seqeunce in each 
# generation, split into herb, crop, plow and spot
function get_best_seq(sol_ob::Tuple, A::Tuple)

	best_ind = best_n_seq(1, sol_ob[2])

	N_gen = length(sol_ob[1])
	N_t = size(sol_ob[1][1])[2]

	herb = Array{Int64, 2}(N_gen, N_t)
	crop = Array{Int64, 2}(N_gen, N_t)
	plow = Array{Int64, 2}(N_gen, N_t)
	spot = Array{Int64, 2}(N_gen, N_t)

	for g in 1:N_gen

		acts = act_seq_2_sub_act(A, sol_ob[1][g][best_ind[g], :])

		herb[g, :] = acts[:, ACT_HERB]
		crop[g, :] = acts[:, ACT_CROP]
		plow[g, :] = acts[:, ACT_PLOW]
		spot[g, :] = acts[:, ACT_SPOT]

	end

	return (herb, crop, plow, spot)

end

# show best single sequence, formatted for compact display using ggplot2
# in R for easier faceting and labelling options (longform data.frame)
function best_seq_2_df(sol_ob::Tuple, A::Tuple)

	best_ind = best_n_seq(1, sol_ob[2])

	N_t = size(sol_ob[1][1])[2]

	acts = act_seq_2_sub_act(A, sol_ob[1][end][best_ind[end], :])

	# make into long form dataframe for plotting 
	return DataFrame(sub_act = @data(vcat(repmat(["herbicide"], N_t), 
			repmat(["crop"], N_t), repmat(["plow"], N_t), 
			repmat(["spot"], N_t))),
		ts = @data(repmat(1:N_t, 4)),
		best_act = @data(vcat(acts[:, ACT_HERB], 
			acts[:, ACT_CROP], acts[:, ACT_PLOW],
			acts[:, ACT_SPOT])))

end

# simulate a given action seqence to get the seed bank size and resistance 
# at time step to then map that to the action seqence 
function sim_act_seq(herb_seq::Array{Int64, 1}, crop_seq::Array{Int64, 1}, 
		spot_seq::Array{Int64, 1}, plow_seq_int::Array{Int64, 1}, 
		pars::DataFrame, low_g::Float64, up_g::Float64, dg::Float64)

	T = length(herb_seq)

	#A = make_action_space()

	g_vals = collect(low_g : dg : up_g)
	len_g = size(g_vals)[1]
	g1_vals = repeat(g_vals, inner = len_g)
	g2_vals = repeat(g_vals, outer = len_g)

	Ng = length(g1_vals)

	# pre-calc the effet of managment actions 
	crop_sur_tup = (1.0, pars[:sur_alt][1], 0.0)
	spot_sur_tup = (0.0, 1.0)

	herb_sur_tup = survial_herb_setup(g1_vals, g2_vals, pars[:p_ex][1], 
		pars[:s0][1], pars[:eff_h1][1], pars[:eff_h2][1], 
		pars[:p_g1h1][1], pars[:p_g1h2][1], pars[:p_g2h1][1], 
		pars[:p_g2h2][1])

	fit_cost = fec_cost_maker(pars[:fr][1], pars[:f0][1], g1_vals, g2_vals)
	
	# generate the indexing keys
	mix_keys = make_index_keys(len_g, len_g)
	
	# define resistance trait values and their co-var
	cov_mat = [pars[:off_sd][1] pars[:off_cv][1];
		   pars[:off_cv][1] pars[:off_sd][1]]
	mix_kernel = offspring_dist_setup(g1_vals, g2_vals, cov_mat, 
		mix_keys)

	# seedbanks T by Ng in size
	SB1 = zeros(T + 1, Ng) 
	SB2 = zeros(T + 1, Ng) 
	# above ground populations 
	ab_pop = zeros(T + 1, Ng)
	ab_pop_spot = zeros(T + 1, Ng)
	# temp holding, rewritten eacbh time step
	mat_pop = zeros(Ng)
	pat_pop = zeros(Ng)
	eff_pop = zeros(Ng)
	par_mix = zeros(length(mix_keys[1]))

	# intial populations 
	int_cov = [pars[:int_sd][1] pars[:int_cv][1]; 
		   pars[:int_cv][1] pars[:int_sd][1]]
	int_mu = [pars[:int_g1][1]; pars[:int_g2][1]]
	int_dist = MvNormal(int_mu, int_cov);
	int_sb1 = pdf(int_dist, transpose(hcat(g1_vals, g2_vals))) * 
		pars[:int_N][1]  
	int_sb2 = zeros(size(g1_vals)[1]) 

	plow_bool = [false, true]
	
	plow_seq = plow_bool[plow_seq_int]

	one_run!(SB1, SB2, ab_pop, ab_pop_spot, mat_pop,
		pat_pop, eff_pop, par_mix, mix_keys, mix_kernel,
		fit_cost, crop_sur_tup, herb_sur_tup, crop_seq, 
		spot_seq, plow_seq, herb_seq, T, int_sb1, int_sb2, 
		pars[:sur_alt][1], pars[:inv_frac][1], pars[:germ_prob][1],
		pars[:seed_sur][1], pars[:fec_max][1], pars[:fec_dd][1], dg)


	return (SB1, SB2)

end

# summary of populatin and managment performance over time 
function get_SB_size(pop_sim::Tuple, dg::Float64)

	return (sum(pop_sim[1], 2) * dg, sum(pop_sim[2], 2) * dg)

end

function get_sur_herb(pop_sim::Tuple{Array{Float64, 2}, Array{Float64, 2}},
	pars::DataFrame, dg::Float64, low_g::Float64, up_g::Float64)

	tot_pop = get_SB_size(pop_sim, dg)[1]

	T = length(tot_pop)

	g_vals = collect(low_g : dg : up_g)
	len_g = size(g_vals)[1]
	g1_vals = repeat(g_vals, inner = len_g)
	g2_vals = repeat(g_vals, outer = len_g)

	herb_sur_tup = survial_herb_setup(g1_vals, g2_vals, pars[:p_ex][1], 
		pars[:s0][1], pars[:eff_h1][1], pars[:eff_h2][1], 
		pars[:p_g1h1][1], pars[:p_g1h2][1], pars[:p_g2h1][1], 
		pars[:p_g2h2][1])

	sur_herb1 = zeros(T)
	sur_herb2 = zeros(T)
	sur_herb12 = zeros(T)

	for i in 1:T

		if tot_pop[i] > 0

			sur_herb1[i] = sum(pop_sim[1][i, :] .* 
					   herb_sur_tup[2]) * dg
			sur_herb2[i] = sum(pop_sim[1][i, :] .* 
					   herb_sur_tup[3]) * dg
			sur_herb12[i] = sum(pop_sim[1][i, :] .*
					    herb_sur_tup[4]) * dg

		end

	end
	
	return DataFrame(herb1 = vcat(sur_herb1 ./ tot_pop...), 
		herb2 = vcat(sur_herb2 ./ tot_pop...), 
		herb12 = vcat(sur_herb12 ./ tot_pop...))

end



