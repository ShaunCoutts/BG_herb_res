# some helper functions to process the output object of the sensitivity analsyis
# package it up as a dataframe for export to R.

using DataFrames

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"
data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out"

# make a summary of a single run for TSR from sense analysis 
function make_sum_sense(sol_list::Array{Any, 1}, T_used::Int64)

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

		herb_seq = best_seq[:, ACT_HERB]
		crop_seq = best_seq[:, ACT_CROP]
		plow_seq_int = best_seq[:, ACT_PLOW]
		spot_seq = best_seq[:, ACT_SPOT]

		N_G = 9 # number of genotypes 
		mix_key = make_TSR_mix_index(N_G)

		# re-build the populaiton at different snapshots
		sim_pop = sim_act_seq(herb_seq, crop_seq, spot_seq, 
			plow_seq_int, sol, mix_key)

		sim_noact = sim_act_seq(repeat([HERB0], inner = T_used),
			repeat([CROP_WW], inner = T_used), 
			repeat([SPOT0], inner = T_used),
			repeat([PLOW0], inner = T_used), 
			sol, mix_key)
		
		ag_pop = get_above_ground(sim_pop, sol, herb_seq, crop_seq, spot_seq, 
			mix_key)

		#make the summaries	
		df_temp[:proWW] = get_pro_WW(crop_seq, T_used)

		df_temp[:proPlow] = get_pro_plow(plow_seq_int, T_used)

		df_temp[:spot_spend] = get_spot_spend(ag_pop[:pre], sol, 
			spot_seq, T_used)

		df_temp[:pro_max_reward] = get_reward_pro_max(ag_pop, sol, 
			sub_acts, T_used)

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

		print("$i\n")

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

# take the simed population and infer above ground population from it
function get_above_ground(sim_pop::Tuple{Array{Float64, 2}, Array{Float64, 2}},
	pars::Dict, herb_seq::Array{Int64, 1}, crop_seq::Array{Int64, 1}, 
	spot_seq::Array{Int64, 1}, mix_key::Tuple)
	
	T = size(sim_pop[1])[1]

	herb_sur = make_herb_sur_dom(mix_key, pars[:s0], pars[:p_ex], 
		pars[:sur_herb])

	ab_pre_spot = sim_pop[1] * pars[:germ_prob]
	ab_post_spot = zeros(size(ab_pre_spot))

	for t in 2:T

		if crop_seq[t - 1] == CROP_FAL
	  
			ab_pre_spot[t, :] .= 0.0
	  
		elseif crop_seq[t - 1] == CROP_WW
	    
			ab_pre_spot[t, :] = ab_pre_spot[t, :] .* 
				herb_sur[herb_seq[t - 1]]
	    
		else
	    
			ab_pre_spot[t, :] = (ab_pre_spot[t, :] .* 
				herb_sur[herb_seq[t - 1]]) * 
				pars[:sur_alt]
	  
		end

		if spot_seq[t - 1] == 0

			ab_post_spot[t, :] = deepcopy(ab_pre_spot[t, :])  

		end

 	end
  
	return Dict(:pre => ab_pre_spot, :post => ab_post_spot)

end

function get_spot_spend(pre_spot_pop::Array{Float64, 2}, pars::Dict, 
	spot_seq::Array{Int64, 1}, T_used::Int64)

	spot_spend = zeros(T_used)

	for t in 1:T_used

		if spot_seq[t] == 1

			spot_spend[t] = pars[:spot_fix] + 
				pars[:spot_var] * sum(pre_spot_pop[t + 1, :]) 

		end

 	end
  
	return sum(spot_spend)

end

function get_undis_reward(ag_pop::Dict{Symbol, Array{Float64, 2}}, 
	pars::Dict{Symbol, Float64}, sub_acts::Array{Int64, 2})

	# total population sizes 
	tot_ab_pre = vcat(sum(ag_pop[:pre], 2)...)	
	tot_ab_pre = tot_ab_pre[2:end]

	tot_ab_post = vcat(sum(ag_pop[:post], 2)...)
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

function get_reward_pro_max(ag_pop::Dict{Symbol, Array{Float64, 2}}, 
	pars::Dict{Symbol, Float64}, sub_acts::Array{Int64, 2}, T_used::Int64)

	# get undiscounted reward from action 
	undis_rewards = get_undis_reward(ag_pop, pars, sub_acts)	

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
function get_resist_profile(sim_pop::Tuple{Array{Float64, 2}, Array{Float64, 2}},
	pars::Dict{Symbol, Float64}, mix_key::Tuple, t::Int64) 

	herb_sur = make_herb_sur_dom(mix_key, pars[:s0], pars[:p_ex], 
		pars[:sur_herb])

	res_H1 = herb_sur[HERB1] .== pars[:s0]
	res_H2 = herb_sur[HERB2] .== pars[:s0]
	res_H12 = herb_sur[HERB12] .== pars[:s0]

	pop_t = deepcopy(sim_pop[1][t, :])

	tot_pop = sum(pop_t)

	return Dict(:res0 => sum(pop_t[!res_H1 & !res_H2]) / tot_pop,
		    :resH1 => sum(pop_t[res_H1]) / tot_pop,
		    :resH2 => sum(pop_t[res_H2]) / tot_pop,
		    :resH12 => sum(pop_t[res_H12]) / tot_pop)

end

function get_SB_t(sim_pop::Tuple{Array{Float64, 2}, Array{Float64, 2}}, 
	t::Int64)

	return sum(sim_pop[1][1:t, :])

end

# get the seed bank size up to t, reative to the seedbank if no action taken
function get_realtive_SB(sim_pop::Tuple{Array{Float64, 2}, Array{Float64, 2}}, 
	sim_noact::Tuple{Array{Float64, 2}, Array{Float64, 2}}, t::Int64)

	return sum(sim_pop[1][2:(t + 1), :]) / sum(sim_noact[1][2:(t + 1), :])

end








