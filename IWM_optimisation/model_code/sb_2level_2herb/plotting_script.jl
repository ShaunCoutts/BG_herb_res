## Plotting script to visulise the population and optimisation outputs
using Plots
using StatPlots
pyplot()
using HDF5, JLD

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
cd(data_loc)
sol = load("test_sol_bad_int.jld")["test_sol"];

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
plt = plot(layout = lm, size = (800, 1200));

# add the seed bank  
plot!(plt, 0:20, [SB[1] SB[2]], guidefont = Plots.font(14), 
	labels = ["seed bank top" "seed bank bottom"], yguide = "amount", 
	tickfont = Plots.font(12), markershape = :circle, linewidth = 2, 
	seriescolor = [plt_greys[2] plt_greys[1]], markerstrokewidth = 0, 
	markersize = 5, xlims = (0, 20.1), legendfont = Plots.font(12),
	subplot = 1);

# add the resistance plot
plot!(plt, 0:20, [resist[:herb1] resist[:herb2]], guidefont = Plots.font(14), 
	labels = ["herb 1" "herb 2"], yguide = "survival exposed",
	tickfont = Plots.font(12), markershape = :circle, linewidth = 2, 
	seriescolor = [col_pal[2] col_pal[3]], markerstrokewidth = 0, 
	markersize = 5, xlims = (0, 20.1), legendfont = Plots.font(12),
	subplot = 2);

# show the reward undiscounted so not confounded with discount rate 
plot!(plt, 1:20, reward_t, xlims = (0, 20.1), guidefont = Plots.font(14), 
	legend = :none, yguide = "reward (£)", linewidth = 2, 
	tickfont = Plots.font(12), markershape = :circle,
	seriescolor = plt_greys[2], markersize = 5,
	markerstrokewidth = 0, subplot = 3);

# add the colmat best action found	
plot_colmat!(plt, act_colmat, subplot = 4)
plot!(plt, subplot = 4); # second call needed to put the generated plot in scope

cd(plot_loc)
savefig("QD_high_int_state.pdf")

##########################################################################
# take the results from the parameter sweep and explore them a bit, make 
# sure everything worked
cd(data_loc)
sol_list = load("sol_sweep_hipop.jld")["sol_sweep"];

# make plots of rewards 
n_best = 10;

# get the number of plots in a square, not very space efficent but fine 
n_pars = length(sol_list);
gd = get_grid_dim(n_pars);

lm = @layout grid(gd[1], gd[2]);

plt = plot(layout = lm, xlabel = "generation", ylabel = "reward",
	   size = (gd[2] * 500, gd[1] * 500));

for i in 1:n_pars

	sol = sol_list[i];

	pars = sol[3];

	# unpack some of the parameters for printing
	int_g2 = round(pars[:int_g2], 3);
	int_g1 = round(pars[:int_g1], 3);
	int_N = pars[:int_N];
	off_cv = pars[:off_cv];
	dr = pars[:dis_rate];
	y0 = pars[:Y0];

	re_gen = reward_gen(sol[2], best_n_seq(n_best, sol[2]));

	plot!(plt, re_gen, label = reshape(repmat([""], n_best), 1, n_best),
	     subplot = i, 
	     title = "int_g1 = $int_g1|int_g2 = $int_g2|int_N = $int_N\noff_cv = $off_cv|dis_rate = $dr|Y0 = $y0 ");

end

cd(plot_loc)
savefig("rewards_over_gen_hipop.pdf")

# show the whole populaiton of actions in the final population 
new_dir = string(plot_loc, "/sol_pops_hipop")
rm(new_dir, recursive = true)
mkdir(new_dir)
cd(new_dir)

n_pars = length(sol_list)

A = make_action_space()

for i in 1:n_pars

	sol = sol_list[i];

	pars = sol[3];

	# unpack some of the parameters for printing
	int_g2 = round(pars[:int_g2], 3);
	int_g1 = round(pars[:int_g1], 3);
	int_N = pars[:int_N];
	off_cv = pars[:off_cv];
	dr = pars[:dis_rate];
	y0 = pars[:Y0];

	best_seq = get_best_seq(sol, A);

	par_title = "int_g1 = $int_g1|int_g2 = $int_g2\nint_N = $int_N|off_cv = $off_cv\ndis_rate = $dr|Y0 = $y0"

	lm = @layout grid(4, 1);

	plt = plot(layout = lm, size = (500, 1000), 
		xlabel = ["" "" "" "time"], ylabel = ["Gen" "Gen" "Gen" "Gen"]);
	heatmap!(plt, best_seq[ACT_HERB], title = string(par_title, "\nHERB"), subplot = 1);
	heatmap!(plt, best_seq[ACT_CROP], title = "CROP", subplot = 2);
	heatmap!(plt, best_seq[ACT_PLOW], title = "PLOW", subplot = 3);
	heatmap!(plt, best_seq[ACT_SPOT], title = "SPOT", subplot = 4);

	savefig(plt, "sol_pop_$i.pdf")

end

# Visulise the solutions 
new_dir = string(plot_loc, "/sol_vis")
rm(new_dir, recursive = true)
mkdir(new_dir)
cd(new_dir)

n_pars = length(sol_list);

A = make_action_space();

col_pal = make_col_pal();
plt_greys = colormap("Grays", 10)[[5, 10]];
low_g = -10.0;
up_g = 10.0;
dg = 1.0;

for i in 1:n_pars

	plt = solution_viz(sol_list[i], A, up_g, low_g, dg)	
	
	savefig(plt, "sol_vis_$i.pdf")

end

# look at the diversity of solutions in last generation 
# plot the sub actions at over the whole population,  make sure there is 
# good diversity for at least a good number of generations 
new_dir = string(plot_loc, "/sol_diversity")
rm(new_dir, recursive = true)
mkdir(new_dir)
cd(new_dir)

n_pars = length(sol_list);

A = make_action_space();

for i in 1:n_pars
	
	sol = sol_list[i];

	sub_acts = actseq_subact(sol[1][end], A);

	lm = @layout [b{0.0h}; grid(2, 2)];

	pars = sol[3];

	# unpack some of the parameters for printing
	int_g2 = round(pars[:int_g2], 3);
	int_g1 = round(pars[:int_g1], 3);
	int_N = pars[:int_N];
	off_cv = pars[:off_cv];
	dr = pars[:dis_rate];
	y0 = pars[:Y0];

	par_title = "int_g1 = $int_g1|int_g2 = $int_g2\nint_N = $int_N|off_cv = $off_cv\ndis_rate = $dr|Y0 = $y0";

	plt = plot(layout = lm, title = [par_title "HERB" "CROP" "PLOW" "SPOT"], 
		xlabel = ["" "" "" "time" "time"], ylabel = ["" "pop of seq" "" "pop of seq" ""]);
	heatmap!(plt, sub_acts[ACT_HERB], subplot = 2);
	heatmap!(plt, sub_acts[ACT_CROP], subplot = 3);
	heatmap!(plt, sub_acts[ACT_PLOW], subplot = 4);
	heatmap!(plt, sub_acts[ACT_SPOT], subplot = 5);

	savefig(plt, "sol_div_$i.pdf")

end















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

