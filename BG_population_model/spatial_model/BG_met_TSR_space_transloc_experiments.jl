# simulation experiments for the TSR and NTSR joint evolution
# to test how TSR can invade into a population with NTSR

using DataFrames

# make a thin wrapper to pass to pmap that unpacks the parameter values 
@everywhere function runner_wrapper(pars::Array{Any, 1})

  #unpack the parameter vlaues 
  run_res = run_scene_trans(pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], pars[7], pars[8],
    pars[9], pars[10], pars[11], pars[12], pars[13], pars[14], pars[15], pars[16], pars[17], pars[18], 
    pars[19], pars[20], pars[21], pars[22], pars[23], pars[24], pars[25], pars[26], pars[27], pars[28], 
    pars[29], pars[30], pars[31], pars[32], pars[33])
    
    return [pars[8:9]; pars[15:16]; pars[18:32]; trans_ex_snapshot(run_res, convert(Float64, pars[2]), pars[34])]
    
end

# need to run each source scenario three times, once in a empty landscape, one in an exposed population
# and one in a herbicide exposed population (can pre calculate all these so don't need to re-caclulate each time 
@everywhere file_loc_func_p = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model" 
@everywhere cd(file_loc_func_p)
@everywhere include("BG_met_TSR_space_exper_funs.jl")
@everywhere include("BG_met_TSR_space_pop_process.jl")
@everywhere include("BG_met_TSR_space_dispersal_functions.jl")
include("spatial_model_plotting_script.jl")
  

###################################################################################################
# script to run the translocation experiments
# parameter values
upper_g = 10;
lower_g = -10;
dg = 0.5;
int_mean_g = 0.0;
int_sd_g = 1.4142;

x_dim = 500; # number of spatial evaluation points, actual landscape size is x_dim * dx
dx = 1.0;

int_num_RR = 0.0;
int_num_Rr = 0.0;
int_num_rr = 10.0; # number of intial seeds at each location for each genoptype, assume only TS susceptible
burnin = 20;
num_iter = 50;

seed_pro_short = 0.4; 
seed_mean_dist_short = 0.5; 
pro_seeds_to_mean_short = 0.4; 
seed_mean_dist_long = 1.5; 
pro_seeds_to_mean_long = 0.4;
scale_pollen = 32.0;
shape_pollen = 3.23; 
offspring_sd = 1.0;
fec_max = 45.0;
fec0 = 10.0;
fec_cost = 0.3;
dd_fec = 0.15;
base_sur = 10.0; 
herb_effect = 20.0; 
g_prot = 2.5; 
pro_exposed = 0.8;
seed_sur = 0.45;
germ_prob = 0.52;
resist_G = ["RR", "Rr"];

# set up the evaluation points for quantitative resistance
g_vals = collect(lower_g : dg : upper_g);   

# set up the plotting color palettes and plotting output locations
sink_col = [HSL(0, 0, 0) HSL(180, 1.0, 0.5) HSL(220, 1.0, 0.5)]
G_col = [HSL(0, 1, 0.7) HSL(280, 1, 0.7) HSL(125, 1, 0.7)] 
output_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output" 

# The set up the set of seeds added after the burnin period
# Four different source populations: low mean_g and low %R, low mean_g and high %R, 
# high mean_g and low %R, high mean_g and high %R

# low mean_g = 0 
inject_g_low = 0.0;
# high mean_g = g such that s = 0.95 for exposed individuals 
inject_g_high = sur_2_g(0.95, herb_effect, base_sur, g_prot); 

inject_sd_g = int_sd_g;

#low and high TSR
int_TSR_low = 0.05;
int_TSR_high = 0.95;

num_inject = 10.0;
inject_locs = [convert(Int64, x_dim / 2)];

# Do the low g and low TSR scenario
g_low_TSR_low = run_scene_trans(g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, int_TSR_low,
  inject_g_low, inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
  resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, g_prot, pro_exposed, 
  seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
  pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd);

plot_scenario(g_low_TSR_low, G_col, sink_col, burnin, output_loc, "g_low_TSR_low_overview.pdf", 1.9)
  
# low g and high TSR source
g_low_TSR_high = run_scene_trans(g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, int_TSR_high,
  inject_g_low, inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
  resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, g_prot, pro_exposed, 
  seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
  pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd);

plot_scenario(g_low_TSR_high, G_col, sink_col, burnin, output_loc, "g_low_TSR_high_overview.pdf", 1.9)

# high g and low TSR source
g_high_TSR_low = run_scene_trans(g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, int_TSR_low,
  inject_g_high, inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
  resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, g_prot, pro_exposed, 
  seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
  pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd);

plot_scenario(g_high_TSR_low, G_col, sink_col, burnin, output_loc, "g_high_TSR_low_overview.pdf", 1.9)

# high g and high TSR source
g_high_TSR_high = run_scene_trans(g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, int_TSR_high,
  inject_g_high, inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
  resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, g_prot, pro_exposed, 
  seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
  pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd);

plot_scenario(g_high_TSR_high, G_col, sink_col, burnin, output_loc, "g_high_TSR_high_overview.pdf", 1.9)


# make some plots that show the effect of a few parameters of interest, try using pmap to speed things 
# up a bit cuase I will have to run these modesl many times

# set a threshold after which we do not bother to calculate spread
threshold = 0.95

# build a parameter list with 1 variable that changes between entries
# start with g_protection
g_prot_var = collect(0 : 1.0 : 4)
par_list_g_pro = []
for i in 1:length(g_prot_var)
  push!(par_list_g_pro, [g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, int_TSR_low,
    inject_g_low, inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
    resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, g_prot_var[i], pro_exposed, 
    seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
    pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd, threshold])
end

# TODO, write the warapper function to unpack this tuple so I can give that to pmap.
# also need that function to return an array or tuple with a couple of summary results

runner_wrapper(par_list_g_pro[1])

out = pmap(runner_wrapper, par_list_g_pro, batch_size = 2)
#TODO: test this pmap function to see if it works





























