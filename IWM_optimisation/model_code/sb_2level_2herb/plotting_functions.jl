using DataFrames
using Colors

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
		pars[:p_g1h1][1], pars[:p_g2h2][1])

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
		pars[:seed_sur][1], pars[:fec_max][1], pars[:fec_dd][1], 
		pars[:sur_spot][1], dg)

	return (SB1, SB2)

end

# summary of populatin and managment performance over time 
function get_SB_size(pop_sim::Tuple, dg::Float64)

	return (sum(pop_sim[1], 2) * dg * dg, 
		sum(pop_sim[2], 2) * dg * dg)

end

function get_sur_herb(pop_sim::Tuple{Array{Float64, 2}, Array{Float64, 2}},
	pars::DataFrame, dg::Float64, low_g::Float64, up_g::Float64)

	tot_pop = get_SB_size(pop_sim, dg)[1]

	T = length(tot_pop)

	g_vals = collect(low_g : dg : up_g)
	len_g = size(g_vals)[1]
	g1_vals = repeat(g_vals, inner = len_g)
	g2_vals = repeat(g_vals, outer = len_g)

	herb_sur_tup = survial_herb_setup(g1_vals, g2_vals, 1.0, 
		pars[:s0][1], pars[:eff_h1][1], pars[:eff_h2][1], 
		pars[:p_g1h1][1], pars[:p_g2h2][1])

	sur_herb1 = zeros(T)
	sur_herb2 = zeros(T)
	sur_herb12 = zeros(T)

	for i in 1:T

		if tot_pop[i] > 0

			sur_herb1[i] = sum(pop_sim[1][i, :] .* 
				herb_sur_tup[2]) * dg * dg
			sur_herb2[i] = sum(pop_sim[1][i, :] .* 
				herb_sur_tup[3]) * dg * dg
			sur_herb12[i] = sum(pop_sim[1][i, :] .*
				herb_sur_tup[4]) * dg * dg

		end

	end
	
	return DataFrame(herb1 = vcat(sur_herb1 ./ tot_pop...), 
		herb2 = vcat(sur_herb2 ./ tot_pop...), 
		herb12 = vcat(sur_herb12 ./ tot_pop...))

end


function get_undis_reward(pop_sim::Tuple{Array{Float64, 2}, Array{Float64, 2}},
		A::Tuple, act_seq::Array{Int64, 1},pars::DataFrame, 
		dg::Float64, low_g::Float64, up_g::Float64)

	T = size(pop_sim[1])[1]

	g_vals = collect(low_g : dg : up_g)
	len_g = size(g_vals)[1]
	g1_vals = repeat(g_vals, inner = len_g)
	g2_vals = repeat(g_vals, outer = len_g)

	herb_sur_tup = survial_herb_setup(g1_vals, g2_vals, pars[:p_ex][1], 
		pars[:s0][1], pars[:eff_h1][1], pars[:eff_h2][1], 
		pars[:p_g1h1][1], pars[:p_g2h2][1])

	ab_pop = pop_sim[1] * pars[:germ_prob][1]
	ab_pop_spot = zeros(size(ab_pop))

	sub_acts = act_seq_2_sub_act(A, act_seq)

	for t in 2:T

		if sub_acts[t - 1, ACT_CROP] == CROP_FAL
	  
			ab_pop[t, :] .= 0.0
	  
		elseif sub_acts[t - 1, ACT_CROP] == CROP_WW
	    
			ab_pop[t, :] = ab_pop[t, :] .* 
				herb_sur_tup[sub_acts[t - 1, ACT_HERB]]
	    
		else
	    
			ab_pop[t, :] = (ab_pop[t, :] .* 
				herb_sur_tup[sub_acts[t - 1, ACT_HERB]]) * 
				sur_crop_alt
	  
		end

		if sub_acts[t - 1, ACT_SPOT] == 1

			ab_pop_spot[t, :] = ab_pop[t, :] * pars[:sur_spot][1]

		end

 	end
  
	tot_ab_pop = vcat(sum(ab_pop, 2)...) * dg * dg
	tot_ab_pop = tot_ab_pop[2:end]

	tot_ab_pop_spot = vcat(sum(ab_pop_spot, 2)...) * dg * dg
	tot_ab_pop_spot = tot_ab_pop_spot[2:end]

	cost_space = make_cost_space(pars[:cost_herb][1], pars[:cost_WW][1],
		pars[:cost_alt][1], pars[:cost_fal][1], pars[:cost_plow][1])
				     
	
	reward = economic_reward(tot_ab_pop_spot, sub_acts[:, ACT_CROP],
			pars[:Y0][1], pars[:Y_slope][1], pars[:Y_alt][1], 
			pars[:rep_pen][1]) - 
		costs(sub_acts, cost_space, pars[:cost_spot][1], tot_ab_pop) 

	return reward

end

##########################################################################
# set out some colour pallets to help consistency 
function make_col_pal()

	bin_act = colormap("Grays", 10)[[1, 10]]

	herb_col = colormap("Blues", 10)[[4, 7, 10]]
	crop_col = colormap("Greens", 10)[[4, 7, 10]]
	plow_col = [bin_act[1], bin_act[2]]
	spot_col = [bin_act[1], bin_act[2]]

	return vcat([bin_act[1], herb_col, crop_col, plow_col, spot_col]...)

end

# re-code the action numbers into sequental numbers that match the above
# pallett
function recode_act(sub_acts::Array{Int64, 2})
	
	recode_act = deepcopy(sub_acts)

	recode_act[:, ACT_CROP] += 4
	recode_act[:, ACT_PLOW] += 7
	recode_act[:, ACT_SPOT] += 10
	
	return recode_act

end

# convert the integer action numbers to RGB colour arrays. also transpose
# and reorder the rows so it is nicer to plot 
function subact_2_colmat(sub_acts::Array{Int64, 2}, 
	col_pal::Array{ColorTypes.RGB{Float64}, 1})

	# recode the actions so they are all on the same pallet
	rc_act = recode_act(sub_acts)

	herb_col = reshape(col_pal[rc_act[:, ACT_HERB]], 1, 
		length(rc_act[:, ACT_HERB]))
	crop_col = reshape(col_pal[rc_act[:, ACT_CROP]], 1, 
		length(rc_act[:, ACT_CROP]))
	plow_col = reshape(col_pal[rc_act[:, ACT_PLOW]], 1, 
		length(rc_act[:, ACT_PLOW]))
	spot_col = reshape(col_pal[rc_act[:, ACT_SPOT]], 1, 
		length(rc_act[:, ACT_SPOT]))

	return vcat([herb_col, crop_col, plow_col, spot_col]...)

end

# make a recatangle function
rect(x, y, w, h) = Shape(x + [0, w, w, 0], y + [0, 0, h, h])

# plot the color map on the right 
function plot_colmat!(plt::Plots.Plot, col_map::Array{ColorTypes.RGB{Float64}, 2};
		    subplot = 1)

	# get the dimentions of the plot
	plt_dims = size(col_map)

	# make the intial plot
	plot!(plt, xlim = (0, plt_dims[2]), ylim = (0, plt_dims[1]), 
	      subplot = subplot, legend = :none)

	# add the rectangles
	for x in 0:(plt_dims[2] - 1)
		for y in 0:(plt_dims[1] - 1)

			plot!(plt, rect(x, plt_dims[1] - y, 1, -1), 
				c = col_map[y + 1, x + 1], 
			      	xlabel = "time step", subplot = subplot,
				xtickfont = Plots.font(9),  
				ytickfont = Plots.font(0))

		end
	end

	ann_lit = colormap("Grays", 10)[7]

	plot!(plt, subplot = subplot,
	      annotations = [(2.4, 3.5, text("Herbicide", ann_lit, 20)),
			(2.4, 2.5, text("Crop       ", ann_lit, 20)),
			(2.4, 1.5, text("Plow       ", ann_lit, 20)),
			(2.4, 0.5, text("Spot       ", ann_lit, 20))])
	
	return nothing

end









