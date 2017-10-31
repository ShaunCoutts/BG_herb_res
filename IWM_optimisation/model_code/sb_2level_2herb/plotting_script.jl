## Plotting script to visulise the population and optimisation outputs
using Plots
using StatPlots
pyplot()


using Colors

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"

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


sb_size = get_SB_size(sim_pop, dg)

surs = get_sur_herb(sim_pop, sol[3], dg, low_g, up_g)





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

