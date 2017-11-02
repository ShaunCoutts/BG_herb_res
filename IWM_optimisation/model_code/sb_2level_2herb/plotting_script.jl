## Plotting script to visulise the population and optimisation outputs
using Plots
using StatPlots
pyplot()

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"
data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out"
plot_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/plot_out"

cd(code_loc)
include("plotting_functions.jl"); 
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

sub_acts = actseq_subact(sol[1][30], A)

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

# simulate the best action seqence to then pull out some metrics of 
# interest from the population 
best_seq = get_best_seq(sol, A);

herb_seq = best_seq[ACT_HERB][end, :];
crop_seq = best_seq[ACT_CROP][end, :];
plow_seq_int = best_seq[ACT_PLOW][end, :];
spot_seq = best_seq[ACT_SPOT][end, :];

sim_pop = sim_act_seq(herb_seq, crop_seq, spot_seq, plow_seq_int, 
		sol[3], low_g, up_g, dg)

# plot the population over time under the best action sequence

# need to get the integer version since some of these functions call 
# internal solver functions 
best_act = sol[1][end][best_n_seq(1, sol[2])[end], :]; 	
best_sub_acts = act_seq_2_sub_act(A, best_act);

reward_t = get_undis_reward(sim_pop, A, best_act, sol[3], dg, low_g, up_g);
resist = get_sur_herb(sim_pop, sol[3], dg, low_g, up_g);
SB = get_SB_size(sim_pop, dg);

# get the action as a col_mat
col_pal = make_col_pal();
act_colmat = subact_2_colmat(best_sub_acts, col_pal);
plt_greys = colormap("Grays", 10)[[5, 10]]

lm = @layout [grid(3, 1)
	     b{0.2h}]
plt = plot(layout = lm)

# add the seed bank  
plot!(plt, 0:20, [SB[1] SB[2]], 
	labels = ["seed bank top" "seed bank bottom"], yguide = "amount", 
	xtickfont = Plots.font(9), markershape = :circle, linewidth = 2, 
	seriescolor = [plt_greys[2] plt_greys[1]], markerstrokewidth = 0, 
	markersize = 5, xlims = (0, 20.1), subplot = 1)

# add the resistance plot
plot!(plt, 0:20, [resist[:herb1] resist[:herb2]], 
	labels = ["herb 1" "herb 2"], yguide = "survival exposed",
	xtickfont = Plots.font(9), markershape = :circle, linewidth = 2, 
	seriescolor = [col_pal[2] col_pal[3]], markerstrokewidth = 0, 
	markersize = 5, xlims = (0, 20.1), subplot = 2)

# show the reward undiscounted so not confounded with discount rate 
plot!(plt, 1:20, reward_t, xlims = (0, 20.1), 
	legend = :none, yguide = "reward (£)", linewidth = 2, 
	xtickfont = Plots.font(9), markershape = :circle,
	seriescolor = plt_greys[2], markersize = 5,
	markerstrokewidth = 0, subplot = 3)

# add the colmat best action found	
plot_colmat!(plt, act_colmat, subplot = 4)
plot!(plt, subplot = 4) # second call needed to put the generated plot in scope

#NOT PERFECT BUT GETTING CLOSE TO WHAT I WANT
# to the other plot. also remove numbers on y plot 3, change colour of lines to 
# be consistent with action make the lines in top plot some version of purple or grey
# make lines thicker. add time step axis label. Also timing is off by a bit


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

