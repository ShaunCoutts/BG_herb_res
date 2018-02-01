# Proccess the sensitvity analysis results to a data frame
using HDF5, JLD

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"
data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out"

cd(code_loc)
include("pop_process_TSR.jl"); 
include("act_seq_handelers.jl")
include("managment_functions.jl"); 
include("sense_helper_functions.jl")
include("plotting_functions_TSR.jl")

cd(data_loc)
sol_sweep = load("sens_runs_obj.jld")["sens_obj"];

# take the summary from the first 20 years to avoid the end
sense_sum = make_sum_sense(sol_sweep, 20);

# write the results to file
cd(data_loc)
writetable("sense_summary.csv", sense_sum);

##########################################################################
# SO I SORT OF RESOVLED THE MYSTERY, THERE WAS A PROBLEM WITH THE ABOVE 
# GROUND SIMULATION SO NEED TO CHANGE THE SUMARY FUNCTION TO USE THE NEW 
# SIMULATION FUNCTION THAT WORKS

# have a look at some of the indexes that have negative rewards indicies 
# to get to the bottom of why
#sol = deepcopy(sol_sweep[1466])

#A = make_action_space();	
#N_G = 9; # number of genotypes 
#mix_key = make_TSR_mix_index(N_G);
#mix_kern = make_TSR_kernel(mix_key);
#
#act_seq = sol[:act_seq]
#
## take only the parameters by removing the action seqence
#delete!(sol, :act_seq)
#sol = convert(Dict{Symbol, Float64}, sol)
#
## get the parameter values and put them in data frame
#df_temp = DataFrame(; sol...)
#
## calculate the measures and add them to df
#best_seq = act_seq_2_sub_act(A, act_seq); 
#
#herb_seq = best_seq[:, ACT_HERB];
#crop_seq = best_seq[:, ACT_CROP];
#plow_seq_int = best_seq[:, ACT_PLOW];
#spot_seq = best_seq[:, ACT_SPOT];
#
#sim_pop = sim_act_seq(herb_seq, crop_seq, spot_seq, 
#	plow_seq_int, sol, mix_kern, mix_key);
#
#ag_pop = get_above_ground(sim_pop, sol, herb_seq, crop_seq, spot_seq, mix_key)
#
## dig into the reward function a bit more
#tot_ab_pre = vcat(sum(ag_pop[:pre], 2)...)	
#tot_ab_pre = tot_ab_pre[2:end]
#
#tot_ab_post = vcat(sum(ag_pop[:post], 2)...)
#tot_ab_post = tot_ab_post[2:end]
#
#pars = deepcopy(sol)
#
#cost_space = make_cost_space(pars[:cost_herb], pars[:cost_WW],
#		pars[:cost_alt], pars[:cost_fal], pars[:cost_plow])
#				     
## calcualte the reward
#reward = economic_reward(tot_ab_post, best_seq[:, ACT_CROP],
#		pars[:Y0], pars[:Y_slope], pars[:Y_alt], 
#		pars[:rep_pen]) - 
#	costs(best_seq, cost_space, pars[:spot_fix], pars[:spot_var], 
#		tot_ab_post) 
#
#
#	pars[:Y0] - (pars[:Y_slope] * sum(ag_pop[:post], 2))
#
#p = merge(pars, Dict(:T => 25, :pop_size => 1000, :num_gen => 100, :mut => 0.03))
#
## try and re-run this single parameter set to see if it gets the same weird answer.
#sol2 = GA_solve_verbose_TSR(convert(Int64, p[:T]), convert(Int64, p[:pop_size]), 
#		convert(Int64, p[:num_gen]), p[:mut], p[:cost_herb], p[:cost_WW], 
#		p[:cost_alt], p[:cost_fal], p[:cost_plow], p[:spot_fix], p[:spot_var], 
#		p[:sur_alt], p[:sur_herb], p[:int_N], p[:int_RR], p[:int_Rr], 
#		p[:int_AA], p[:int_Aa], p[:inv_frac], p[:germ_prob], 
#		p[:seed_sur], p[:fec_max], p[:fec_dd], p[:sur_spot], 
#		p[:dis_rate], p[:Y0], p[:Y_slope],p[:Y_alt], p[:p_ex], 
#		p[:s0], p[:rep_pen]);
#	
#
#best_ind = best_n_seq(1, sol2[2])
#out = merge(sol2[3], Dict(:act_seq => sol2[1][end][best_ind[end], :]))
#
## reward from the solver 
#sol2[2][end][best_ind[end], :] # 3773.5, not very large but not negative dig into this a bit more
#
## run the new solver solution through the population reconstruction to see where they diverge
#act_seq2 = out[:act_seq]
#
#delete!(out, :act_seq)
#out = convert(Dict{Symbol, Float64}, out)
#
#best_seq2 = act_seq_2_sub_act(A, act_seq2); 
#
#herb_seq2 = best_seq2[:, ACT_HERB];
#crop_seq2 = best_seq2[:, ACT_CROP];
#plow_seq_int2 = best_seq2[:, ACT_PLOW];
#spot_seq2 = best_seq2[:, ACT_SPOT];
#
#sim_pop2 = sim_act_seq(herb_seq2, crop_seq2, spot_seq2, 
#	plow_seq_int2, out, mix_kern, mix_key);
#
## take a look at both the simulated pop from the original solve and the new one
#hcat(sum(sim_pop[1], 2), sum(sim_pop2[1], 2)) # basically identical
#
## now look at the seed banks the solver was recording 
#sol_SB1 = sol2[4][:SB1_fin][best_ind[end]];
#
#hcat(sum(sim_pop[1], 2), sum(sim_pop2[1], 2), sum(sol_SB1, 2)) # all identical, so the simulation of the action sequence is working  
#
## look at the way the reward is created instead 
#ag_pop2 = get_above_ground(sim_pop2, out, herb_seq2, crop_seq2, spot_seq2, mix_key)
#
#tot_ab_pre2 = vcat(sum(ag_pop2[:pre], 2)...)	
#tot_ab_pre2 = tot_ab_pre2[2:end]
#
#tot_ab_post2 = vcat(sum(ag_pop2[:post], 2)...)
#tot_ab_post2 = tot_ab_post2[2:end]
#
#cost_space2 = make_cost_space(out[:cost_herb], out[:cost_WW],
#		out[:cost_alt], out[:cost_fal], out[:cost_plow])
#				     
## calcualte the reward
#reward2 = economic_reward(tot_ab_post2, best_seq2[:, ACT_CROP],
#		out[:Y0], out[:Y_slope], out[:Y_alt], 
#		out[:rep_pen]) - 
#	costs(best_seq2, cost_space2, out[:spot_fix], out[:spot_var], 
#		tot_ab_post2) 
#
#comp_reward = hcat(out[:Y0] - (out[:Y_slope] * sum(ag_pop2[:post], 2)), 
#     pars[:Y0] - (pars[:Y_slope] * sum(ag_pop[:post], 2)))
#sum(comp_reward, 1) # this gets no where near 3000 that we see in the reward, from the solver 
#
## try using the exact same reward function used in solver with these inputs 
#dis_rates = make_dis_rate(p[:T], out[:dis_rate])
#reward_TSR(sim_pop2[:ag_pre], sim_pop2[:ag_post], dis_rates, out[:Y0], out[:Y_slope], 
#	out[:Y_alt], out[:rep_pen], best_seq2 , cost_space2, 
#	out[:spot_fix], out[:spot_var]) 
#
## -2.3k so very different from the +3k seen in the solver, seedbank is identical so must be a 
## problem in the conversion to the aboveground population somewhere
## compaire the simulated version with the version seen by the solver
#tracked_ad_pre = sol2[4][:ab_pre_fin][best_ind[end]];
#hcat(sum(sim_pop2[:ag_pre], 2), sum(tracked_ad_pre, 2))
#

#these are very different, need to dig into the get_above_ground() function 

# Still no idea, I might need to redo the solver to also produce the population 
# the solver sees.
