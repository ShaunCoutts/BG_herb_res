## Plotting script to visulise the population and optimisation outputs
using Plots 
pyplot()

using Colors
using DataFrames

data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out/"

plot_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/plot_out"

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
function sim_act_seq(act_seq::Array{Int64, 1}, A::Tuple)

	# TODO: MAKE THE SIMULATION, ESSENTIALLY JUST CALL 
	# ONRERUN! AGAIN, WHICH WILL NEED ALL THE PARAMETER 
	# VALUES USED AGAIN. SO MAKE A CALLING WRAPPER THAT RUNS THE
	# MODEL AND PACKAGES UP THE RESULTS WITH THE PARAMETERS


end

##########################################################################
# plots of different metrics and behaviour. Ultimatley ggplot2 will be
# more flexible. So make functions that will generate longform dataframes
# which is then passed to a R plotting function


n_plt = 10

re_gen = reward_gen(sol[2], best_n_seq(n_plt, sol[2]))

plot(re_gen, label = reshape(repmat([""], n_plt), 1, n_plt), 
    xlabel = "generation", ylabel = "reward")


# plot the sub actions at over the whole population,  make sure there is 
# good diversity for at least a good number of generations 
A = make_action_space()

sub_acts = actseq_subact(sol[1][20], A)

lm = @layout grid(2, 2)

plt = plot(layout = lm, title = ["HERB" "CROP" "PLOW" "SPOT"], 
	   xlabel = ["" "" "time" "time"], ylabel = ["pop of seq" "" "pop of seq" ""])
heatmap!(plt, sub_acts[ACT_HERB], subplot = 1)
heatmap!(plt, sub_acts[ACT_CROP], subplot = 2)
heatmap!(plt, sub_acts[ACT_PLOW], subplot = 3)
heatmap!(plt, sub_acts[ACT_SPOT], subplot = 4)

# plot of the best sequences over time to visulise how they change 

best_seq = get_best_seq(sol, A)

lm = @layout grid(4, 1)

plt = plot(layout = lm, title = ["HERB" "CROP" "PLOW" "SPOT"], 
	   xlabel = ["" "" "" "time"], ylabel = ["Gen" "Gen" "Gen" "Gen"])
heatmap!(plt, best_seq[ACT_HERB], subplot = 1)
heatmap!(plt, best_seq[ACT_CROP], subplot = 2)
heatmap!(plt, best_seq[ACT_PLOW], subplot = 3)
heatmap!(plt, best_seq[ACT_SPOT], subplot = 4)


df = best_seq_2_df(sol, A)
# save for shipping to R
cd(data_loc)
writetable("test_GA.csv", df)



# make some colours 
crop_col = colormap("Blues", 4)

#heatmap to plot population over 2d
heatmap(matrix)


# to animate a plot over time
f0s = convert(Array{Float64}, 1:10);
fr = 0.5
an_test = zeros(10, size(g1_vals)[1])

for i = 1:10 

	an_test[i, :] = 1 ./ fec_cost_maker(fr, f0s[i], g1_vals, g2_vals);

end

anim = @animate for i = 1:10
  
	plot(an_test[i,:])
 
end

gif(anim, "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/test.gif", 
  fps = 1)

