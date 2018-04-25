# some helper functions to process the output object of the sensitivity analsyis
# package it up as a dataframe for export to R.

using DataFrames

# take a set of parameters and a action seqence and simulate the population
# then put in data frame to join with the acrtion seq one
function make_pop_df(sol_list::Array{Dict{Symbol, Any}, 1}, T_used::Int64)

	df_list = []
	np = length(sol_list)

	A = make_action_space()

	N_G = 9 # number of genotypes
	mix_key = make_TSR_mix_index(N_G)

	mix_kern = make_TSR_kernel(mix_key)

	for i in 1:np

		sol = deepcopy(sol_list[i])

		# get the action sequence
		act_seq = sol[:act_seq]

		# take only the parameters by removing the action seqence
		delete!(sol, :act_seq)
		sol = convert(Dict{Symbol, Float64}, sol)

		# get the parameter values and put them in data frame
		df_temp = deepcopy(DataFrame(; sol...));
		dfr = deepcopy(df_temp);
		# replicate these parameter combinatins as to make it the right length
		for i in 1:T_used

			append!(dfr, df_temp);

		end

		# calculate the measures and add them to df
		best_seq = act_seq_2_sub_act(A, act_seq)

		herb_seq = best_seq[1:T_used, ACT_HERB]
		crop_seq = best_seq[1:T_used, ACT_CROP]
		plow_seq_int = best_seq[1:T_used, ACT_PLOW]
		spot_seq = best_seq[1:T_used, ACT_SPOT]

		# re-build the populaiton at different snapshots
		sim_pop = sim_act_seq(herb_seq, crop_seq, spot_seq,
			plow_seq_int, sol, mix_kern, mix_key)

		dfr[:ts] = 0:T_used
		dfr[:herb] = vcat([0], herb_seq)
		dfr[:crop] = vcat([0], crop_seq)
		dfr[:plow] = vcat([0], plow_seq_int)
		dfr[:spot] = vcat([0], spot_seq)
		dfr[:SB1] = sum(sim_pop[:SB1], 2)[1:end]
		dfr[:SB2] = sum(sim_pop[:SB2], 2)[1:end]
		dfr[:ag_post] = sum(sim_pop[:ag_post], 2)[1:end]
		dfr[:ag_pre] = sum(sim_pop[:ag_pre], 2)[1:end]
		dfr[:RRAA] = sim_pop[:SB1][1:T_used + 1, 1]
		dfr[:RRAa] = sim_pop[:SB1][1:T_used + 1, 2]
		dfr[:RRaa] = sim_pop[:SB1][1:T_used + 1, 3]
		dfr[:RrAA] = sim_pop[:SB1][1:T_used + 1, 4]
		dfr[:RrAa] = sim_pop[:SB1][1:T_used + 1, 5]
		dfr[:Rraa] = sim_pop[:SB1][1:T_used + 1, 6]
		dfr[:rrAA] = sim_pop[:SB1][1:T_used + 1, 7]
		dfr[:rrAa] = sim_pop[:SB1][1:T_used + 1, 8]
		dfr[:rraa] = sim_pop[:SB1][1:T_used + 1, 9]

		push!(df_list, dfr)

	end

	return vcat(df_list...)

end

# make a summary of a single run for TSR from sense analysis
function make_sum_sense(sol_list::Array{Dict{Symbol, Any}, 1}, T_used::Int64)

	df_list = []
	np = length(sol_list)

	A = make_action_space()

	N_G = 9 # number of genotypes
	mix_key = make_TSR_mix_index(N_G)

	mix_kern = make_TSR_kernel(mix_key)

	for i in 1:np

		sol = deepcopy(sol_list[i])

		# get the action sequence
		act_seq = sol[:act_seq]

		# take only the parameters by removing the action seqence
		delete!(sol, :act_seq)
		sol = convert(Dict{Symbol, Float64}, sol)

		# get the parameter values and put them in data frame
		df_temp = DataFrame(; sol...)

		# calculate the measures and add them to df
		best_seq = act_seq_2_sub_act(A, act_seq)

		herb_seq = best_seq[:, ACT_HERB]
		crop_seq = best_seq[:, ACT_CROP]
		plow_seq_int = best_seq[:, ACT_PLOW]
		spot_seq = best_seq[:, ACT_SPOT]

		# re-build the populaiton at different snapshots
		sim_pop = sim_act_seq(herb_seq, crop_seq, spot_seq,
			plow_seq_int, sol, mix_kern, mix_key)

		sim_noact = sim_act_seq(repeat([HERB0], inner = T_used),
			repeat([CROP_WW], inner = T_used),
			repeat([SPOT0], inner = T_used),
			repeat([PLOW0], inner = T_used),
			sol, mix_kern, mix_key)

		#make the summaries
		df_temp[:proWW] = get_pro_WW(crop_seq, T_used)

		df_temp[:proPlow] = get_pro_plow(plow_seq_int, T_used)

		df_temp[:spot_spend] = get_spot_spend(sim_pop[:ag_pre], sol,
			spot_seq, T_used)

		df_temp[:pro_max_reward] = get_reward_pro_max(sim_pop, sol,
			best_seq, T_used)

		df_temp[:herb_apps] = get_num_herb(herb_seq, T_used)

		res_profile = get_resist_profile(sim_pop, sol, mix_key, T_used)

		df_temp[:fin_sucep] = res_profile[:res0]
		df_temp[:fin_resH1] = res_profile[:resH1]
		df_temp[:fin_resH2] = res_profile[:resH2]
		df_temp[:fin_resH12] = res_profile[:resH12]

		df_temp[:SB5] = get_SB_t(sim_pop, 5)
		df_temp[:SB10] = get_SB_t(sim_pop, 10)

		df_temp[:relSB5] = get_realtive_SB(sim_pop, sim_noact, 5)
		df_temp[:relSB10] = get_realtive_SB(sim_pop, sim_noact, 10)

		push!(df_list, df_temp)

		if i % 500 == 0

			print("$i\n")

		end

	end

	return vcat(df_list...)

end

# a set of functions to deal with the out out of TSR sumulation and sense
# stuff
function get_pro_WW(crop_seq::Array{Int64, 1}, T_used::Int64)

	return sum(crop_seq[1:T_used] .== CROP_WW) / T_used

end

function get_pro_plow(plow_seq::Array{Int64, 1}, T_used::Int64)

	return sum(plow_seq[1:T_used] .- 1) / T_used

end

function get_spot_spend(pre_spot_pop::Array{Float64, 2}, pars::Dict,
	spot_seq::Array{Int64, 1}, T_used::Int64)

	spot_spend = zeros(T_used)

	for t in 2:T_used

		if spot_seq[t] == 1

			spot_spend[t] = pars[:spot_fix] +
				pars[:spot_var] * sum(pre_spot_pop[t + 1, :])

		end

 	end

	return sum(spot_spend)

end

function get_undis_reward(sim_pop::Dict{Symbol, Array{Float64, 2}},
	pars::Dict{Symbol, Float64}, sub_acts::Array{Int64, 2})

	# total population sizes
	tot_ab_pre = vcat(sum(sim_pop[:ag_pre], 2)...)
	tot_ab_pre = tot_ab_pre[2:end]

	tot_ab_post = vcat(sum(sim_pop[:ag_post], 2)...)
	tot_ab_post = tot_ab_post[2:end]

	# generate the costs for each action
	cost_space = make_cost_space(pars[:cost_herb], pars[:cost_WW],
		pars[:cost_alt], pars[:cost_fal], pars[:cost_plow])

	# calcualte the reward
	reward = economic_reward(tot_ab_post, sub_acts[:, ACT_CROP],
			pars[:Y0], pars[:Y_slope], pars[:Y_alt],
			pars[:rep_pen]) -
		costs(sub_acts, cost_space, pars[:spot_fix], pars[:spot_var],
			tot_ab_post)

	return reward

end

function get_reward_pro_max(sim_pop::Dict{Symbol, Array{Float64, 2}},
	pars::Dict{Symbol, Float64}, sub_acts::Array{Int64, 2}, T_used::Int64)

	# get undiscounted reward from action
	undis_rewards = get_undis_reward(sim_pop, pars, sub_acts)

	# get max possible reward assuming no black grass
	# note a minor complication here is if the crop repeate i
	# penelty makes cycling the best roation even without Black grass
	WW_reward = sum(vcat(pars[:Y0] - pars[:cost_WW],
		repeat([(pars[:Y0] * pars[:rep_pen]) - pars[:cost_WW]],
			inner = T_used - 1)))
	ALT_reward = sum(vcat(pars[:Y_alt] - pars[:cost_alt],
		repeat([(pars[:Y_alt] * pars[:rep_pen]) - pars[:cost_alt]],
			inner = T_used - 1)))
	best_crop = max(pars[:Y0] - pars[:cost_WW], pars[:Y_alt] - pars[:cost_alt])
	other_crop = min(pars[:Y0] - pars[:cost_WW], pars[:Y_alt] - pars[:cost_alt])
	CYC_reward = repeat([best_crop, other_crop], outer = T_used)
	CYC_reward = sum(CYC_reward[1:T_used])

	best_pos = max(WW_reward, ALT_reward, CYC_reward)

	return sum(undis_rewards[1:T_used]) / best_pos

end

function get_num_herb(herb_seq::Array{Int64, 1}, T_used::Int64)

	herb_used = 0

	for t in 1:T_used

		if herb_seq[t] == HERB1 || herb_seq[t] == HERB2

			herb_used += 1

		end

		if herb_seq[t] == HERB12

			herb_used += 2

		end

	end

	return herb_used

end

# get a summary of the resistance in a population after t years
function get_resist_profile(sim_pop::Dict{Symbol, Array{Float64, 2}},
	pars::Dict{Symbol, Float64}, mix_key::Tuple, t::Int64)

	herb_sur = make_herb_sur_dom(mix_key, pars[:s0], pars[:p_ex],
		pars[:sur_herb])

	res_H1 = herb_sur[HERB1] .== pars[:s0]
	res_H2 = herb_sur[HERB2] .== pars[:s0]
	res_H12 = herb_sur[HERB12] .== pars[:s0]

	pop_t = deepcopy(sim_pop[:SB1][t, :])

	tot_pop = sum(pop_t)

	return Dict(:res0 => sum(pop_t[.!res_H1 .& .!res_H2]) / tot_pop,
		    :resH1 => sum(pop_t[res_H1]) / tot_pop,
		    :resH2 => sum(pop_t[res_H2]) / tot_pop,
		    :resH12 => sum(pop_t[res_H12]) / tot_pop)

end

function get_SB_t(sim_pop::Dict{Symbol, Array{Float64, 2}},
	t::Int64)

	return sum(sim_pop[:SB1][1:t, :])

end

# get the seed bank size up to t, reative to the seedbank if no action taken
function get_realtive_SB(sim_pop::Dict{Symbol, Array{Float64, 2}},
	sim_noact::Dict{Symbol, Array{Float64, 2}}, t::Int64)

	return sum(sim_pop[:SB1][2:(t + 1), :]) / sum(sim_noact[:SB1][2:(t + 1), :])

end

# make a slightly different summary that returns the time series of sub-actions
function make_ts_sense(sol_list::Array{Dict{Symbol, Any}, 1}, T_used::Int64)

	df_list = []
	np = length(sol_list)

	A = make_action_space()

	for i in 1:np

		sol = deepcopy(sol_list[i])

		# get the action sequence
		act_seq = sol[:act_seq]

		# take only the parameters by removing the action seqence
		delete!(sol, :act_seq)
		sol = convert(Dict{Symbol, Float64}, sol)

		# get the parameter values and put them in data frame
		df_temp = DataFrame(; sol...)

		# calculate the measures and add them to df
		best_seq = act_seq_2_sub_act(A, act_seq)

		herb_seq = best_seq[1:T_used, ACT_HERB]
		crop_seq = best_seq[1:T_used, ACT_CROP]
		plow_seq_int = best_seq[1:T_used, ACT_PLOW]
		spot_seq = best_seq[1:T_used, ACT_SPOT]

		# put that action sequence into a time series matrix for LCSS package
		# that is easy to parse
		herb_string = join(["$(a) " for a in herb_seq])
		crop_string = join(["$(a) " for a in crop_seq])
		plow_string = join(["$(a) " for a in plow_seq_int])
		spot_string = join(["$(a) " for a in spot_seq])


		df_temp[:sub_acts] = herb_string * crop_string * plow_string * spot_string
		df_temp[:T_used] = T_used

		push!(df_list, df_temp)

		if i % 5000 == 0

			print("$i\n")

		end

	end

	return vcat(df_list...)

end
