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
    
    final_pop = trans_ex_snapshot(run_res, convert(Float64, pars[2]), pars[34], pars[6])
    pars_block = vcat(fill(transpose([pars[8:9]; pars[15:16]; pars[18:32]]), 3) ...)
    return hcat(pars_block, final_pop) 
    
end

# need to run each source scenario three times, once in a empty landscape, one in an exposed population
# and one in a herbicide exposed population (can pre calculate all these so don't need to re-caclulate each time 
@everywhere file_loc_func_p = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model" 
@everywhere cd(file_loc_func_p)
@everywhere include("BG_met_TSR_space_exper_funs.jl")
@everywhere include("BG_met_TSR_space_runners.jl")
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

x_dim = 250; # number of spatial evaluation points, actual landscape size is x_dim * dx
dx = 1.0;

int_num_RR = 0.0;
int_num_Rr = 0.0;
int_num_rr = 10.0; # number of intial seeds at each location for each genoptype, assume only TS susceptible
burnin = 20;
num_iter = 100;

seed_pro_short = 0.48; 
seed_mean_dist_short = 0.58; 
pro_seeds_to_mean_short = 0.44; 
seed_mean_dist_long = 1.65; 
pro_seeds_to_mean_long = 0.39;
scale_pollen = 32.0;
shape_pollen = 3.32; 
offspring_sd = 1.0;
fec_max = 60.0;
fec0 = 4.0;
fec_cost = 0.45;
dd_fec = 0.15;
base_sur = 10.0; 
herb_effect = 14.0; 
g_prot = 1.5; 
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

#### g_prot plot ###
# set a threshold after which we do not bother to calculate spread
threshold = 0.95

# build a parameter list with 1 variable that changes between entries
# start with g_protection
g_prot_var = collect(0 : 0.1 : 3)
g_pro_rep = repmat(g_prot_var, 4)
#use the function below t make a matrix with all the combinations I want
TSR_reps = repeat([int_TSR_low, int_TSR_high], inner = length(g_prot_var) * 2)
g_reps = repeat(repeat([inject_g_low, inject_g_high], inner = length(g_prot_var)), outer = 2)
inj_scens = hcat([g_pro_rep, TSR_reps, g_reps] ...)

par_list_g_pro = []
for i in 1:size(inj_scens)[1]
  push!(par_list_g_pro, [g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, inj_scens[i, 2],
    inj_scens[i, 3], inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
    resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, inj_scens[i, 1], pro_exposed, 
    seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
    pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd, threshold])
end

out = pmap(runner_wrapper, par_list_g_pro, batch_size = 55)

# Turn the output into a dataframe, which is a bit easier to work with, with all the different 
# parameters changing and metrics
mat_out = vcat(out ...)
df_g_pro = DataFrame(mat_out)
names!(df_g_pro, [:inj_TSR, :inj_g, :seed_sur, :germ_prob, :fec_max, :dd_fec, :fec0, :fec_cost,
  :s0, :herb_effect, :g_pro, :pro_expo, :P_s_seed, :mds_seed, :P_mds_seed, :mdl_seed, :P_mdl_seed, 
  :scale_pollen, :shape_pollen, :scen, 
  :pro_R, :sur_g, :spr_TSR, :spr_RR, :spr_Rr, :spr_rr])
  
# This dataframe takes quiet a bit to produce, save results so I don't need to keep re-runing the model_output
cd(output_loc)
# writetable("g_pro_prerun.csv", df_g_pro, header = true)
df_g_pro = readtable("g_pro_prerun.csv", header = true)

#plot these scenarios
pop_res_4_scen(df_g_pro, :g_pro, [:pro_R, :sur_g], sink_col, 
  output_loc, "effect_g_pro_resist.pdf", 1.9, convert(String, L"$\rho$"), [convert(String, L"$\%R_{100}$"), convert(String, L"$Sur_{rr}$")], int_TSR_low, 
  inject_g_low, int_TSR_high, inject_g_high)

  
### TSR in injected population
TSR_var = collect(0 : 0.025 : 1)  
TSR_rep = repmat(TSR_var, 2)
g_reps = repeat([inject_g_low, inject_g_high], inner = length(TSR_var))
inj_scens = hcat([TSR_rep, g_reps] ...)
 
par_list_TSR_inj = []
for i in 1:size(inj_scens)[1]
  push!(par_list_TSR_inj, [g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, inj_scens[i, 1],
    inj_scens[i, 2], inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
    resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, g_prot, pro_exposed, 
    seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
    pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd, threshold])
end

out = pmap(runner_wrapper, par_list_TSR_inj, batch_size = 28)
 
mat_out = vcat(out ...)
df_TSR_inj = DataFrame(mat_out)
names!(df_TSR_inj, [:inj_TSR, :inj_g, :seed_sur, :germ_prob, :fec_max, :dd_fec, :fec0, :fec_cost,
  :s0, :herb_effect, :g_pro, :pro_expo, :P_s_seed, :mds_seed, :P_mds_seed, :mdl_seed, :P_mdl_seed, 
  :scale_pollen, :shape_pollen, :scen, 
  :pro_R, :sur_g, :spr_TSR, :spr_RR, :spr_Rr, :spr_rr])
  
# This dataframe takes quiet a bit to produce, save results so I don't need to keep re-runing the model_output
cd(output_loc)
#writetable("TSR_inj_prerun.csv", df_TSR_inj, header = true)
df_TSR_inj = readtable("TSR_inj_prerun.csv", header = true)

TSR_inj_var_plot(df_TSR_inj, :inj_TSR, [:pro_R, :sur_g], sink_col, 
  output_loc, "TSR_inj_final.pdf", 1.9, "%RR in transported seeds", ["%R", "survival rr"], inject_g_low, 
  inject_g_high)
  
# Do a big run to see the effect of herbicide effectivness and protective effect 
# Computationally this will be pretty expensive for all the source scenarions
# testing each parameter at 6 levels is a nice trade-off between run time and grain size
# If I want an even effect of herbicide on the probability scale need to work out the effect 
# on the -(logit scale - base sur)
herb_eff_var = [19.21024, 14.587071, 13.888764, 13.474725, 13.177533, 12.944439]
g_pro_var = collect(1.0 : (3.0 - 1.0) / 5 : 3.0)
sur_var = hcat([repeat(herb_eff_var, inner = length(g_pro_var)), repeat(g_pro_var, outer = length(herb_eff_var))] ...)

TSR_reps = repeat([int_TSR_low, int_TSR_high], inner = size(sur_var)[1] * 2)
g_reps = repeat(repeat([inject_g_low, inject_g_high], inner = size(sur_var)[1]), outer = 2)
inj_scens = hcat([repmat(sur_var, 4), TSR_reps, g_reps] ...)

par_list_herb_eff = []
for i in 1:size(inj_scens)[1]
  push!(par_list_herb_eff, [g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, inj_scens[i, 3],
    inj_scens[i, 4], inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
    resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, inj_scens[i, 1], inj_scens[i, 2], pro_exposed, 
    seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
    pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd, threshold])
end

out = pmap(runner_wrapper, par_list_herb_eff, batch_size = 48)
 
mat_out = vcat(out ...)
df_herb_eff = DataFrame(mat_out)
names!(df_herb_eff, [:inj_TSR, :inj_g, :seed_sur, :germ_prob, :fec_max, :dd_fec, :fec0, :fec_cost,
  :s0, :herb_effect, :g_pro, :pro_expo, :P_s_seed, :mds_seed, :P_mds_seed, :mdl_seed, :P_mdl_seed, 
  :scale_pollen, :shape_pollen, :scen, 
  :pro_R, :sur_g, :spr_TSR, :spr_RR, :spr_Rr, :spr_rr])
  
# This dataframe takes quiet a bit to produce, save results so I don't need to keep re-runing the model_output
cd(output_loc)
# writetable("herb_eff_prerun.csv", df_herb_eff, header = true)
df_herb_eff = readtable("herb_eff_prerun.csv", header = true)

# conver the herb_effect to probability in a new coloumn
df_herb_eff[:prob_herb_eff] = logit_sur_2_prob(df_herb_eff[:herb_effect], df_herb_eff[:s0][1]);

# start with the source pop low g, low TSR
df_glow_TSRlow = df_herb_eff[(df_herb_eff[:inj_g] .== inject_g_low) & (df_herb_eff[:inj_TSR] .== int_TSR_low), :]

colormats_2D(df_glow_TSRlow, :g_pro, :prob_herb_eff, :pro_R, output_loc, 
  "herb_eff_g_pro_lowg_lowTSR.pdf", 1.9, "protective effect g (logits)", "survival suceptable", 
  "%R", "quant. res. low, TSR low")

# source pop low g, high TSR 
df_glow_TSRhigh = df_herb_eff[(df_herb_eff[:inj_g] .== inject_g_low) & (df_herb_eff[:inj_TSR] .== int_TSR_high), :]

colormats_2D(df_glow_TSRhigh, :g_pro, :prob_herb_eff, :pro_R, output_loc, 
  "herb_eff_g_pro_lowg_highTSR.pdf", 1.9, "protective effect g (logits)", "survival suceptable", 
  "%R", "quant. res. low, TSR high")
 
# source pop high g , low TSR 
df_ghigh_TSRlow = df_herb_eff[(df_herb_eff[:inj_g] .== inject_g_high) & (df_herb_eff[:inj_TSR] .== int_TSR_low), :]

colormats_2D(df_ghigh_TSRlow, :g_pro, :prob_herb_eff, :pro_R, output_loc, 
  "herb_eff_g_pro_highg_lowTSR.pdf", 1.9, "protective effect g (logits)", "survival suceptable", 
  "%R", "quant. res. high, TSR low")
   
# source pop high g , high TSR 
df_ghigh_TSRhigh = df_herb_eff[(df_herb_eff[:inj_g] .== inject_g_high) & (df_herb_eff[:inj_TSR] .== int_TSR_high), :]

colormats_2D(df_ghigh_TSRhigh, :g_pro, :prob_herb_eff, :pro_R, output_loc, 
  "herb_eff_g_pro_highg_highTSR.pdf", 1.9, "protective effect g (logits)", "survival suceptable", 
  "%R", "quant. res. high, TSR low")
  
### now for survival 
# start with the source pop low g, low TSR
colormats_2D(df_glow_TSRlow, :g_pro, :prob_herb_eff, :sur_g, output_loc, 
  "herb_eff_g_pro_lowg_lowTSR_surg.pdf", 1.9, "protective effect g (logits)", "survival suceptable", 
  "sur_rr", "quant. res. low, TSR low")

# source pop low g, high TSR 
colormats_2D(df_glow_TSRhigh, :g_pro, :prob_herb_eff, :sur_g, output_loc, 
  "herb_eff_g_pro_lowg_highTSR_surg.pdf", 1.9, "protective effect g (logits)", "survival suceptable", 
  "sur_rr", "quant. res. low, TSR high")
 
# source pop high g , low TSR 
colormats_2D(df_ghigh_TSRlow, :g_pro, :prob_herb_eff, :sur_g, output_loc, 
  "herb_eff_g_pro_highg_lowTSR_surg.pdf", 1.9, "protective effect g (logits)", "survival suceptable", 
  "sur_rr", "quant. res. high, TSR low")
   
# source pop high g , high TSR 
colormats_2D(df_ghigh_TSRhigh, :g_pro, :prob_herb_eff, :sur_g, output_loc, 
  "herb_eff_g_pro_highg_highTSR_surg.pdf", 1.9, "protective effect g (logits)", "survival suceptable", 
  "sur_rr", "quant. res. high, TSR low")
  
#######€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€ 
### Different plot for the spatial version of the non-spatial  

# make a thin wrapper to pass to pmap that unpacks the parameter values 
@everywhere function space_time_wrapper(pars::Array{Any, 1})

  #unpack the parameter vlaues 
  run_res = run_scene_trans_space_exp(pars[1], pars[2], pars[3], pars[4], pars[5], pars[6], pars[7], pars[8],
    pars[9], pars[10], pars[11], pars[12], pars[13], pars[14], pars[15], pars[16], pars[17], pars[18], 
    pars[19], pars[20], pars[21], pars[22], pars[23], pars[24], pars[25], pars[26], pars[27], pars[28], 
    pars[29], pars[30], pars[31], pars[32], pars[33])
   
    n_met = size(run_res[1])[1]
    res_block = hcat(repeat(["empty", "naive", "exposed"], inner = n_met), vcat(run_res[1], run_res[2], run_res[3]))
    pars_block = vcat(fill(transpose([pars[8:9]; pars[15:16]; pars[18:33]]), size(res_block)[1]) ...)
    return hcat(pars_block, res_block) 
    
end

upper_g = 20;
lower_g = -20;
dg = 0.5;
int_mean_g = 0.0;
int_sd_g = 1.4142;

x_dim = 200; # number of spatial evaluation points, actual landscape size is x_dim * dx
dx = 1.0;

int_num_rr = 10.0; # number of intial seeds at each location for each genoptype, assume only TS susceptible
burnin = 20;
num_iter = 100;

seed_pro_short = 0.48; 
seed_mean_dist_short = 0.58; 
pro_seeds_to_mean_short = 0.44; 
seed_mean_dist_long = 1.65; 
pro_seeds_to_mean_long = 0.39;
scale_pollen = 32.0;
shape_pollen = 3.32; 
offspring_sd = 1.0;
fec_max = 60.0;
fec0 = 4.0;
fec_cost = 0.45;
dd_fec = 0.15;
base_sur = 10.0; 
herb_effect = 16.0; 
g_prot = 1.5; 
pro_exposed = 0.8;
seed_sur = 0.45;
germ_prob = 0.52;
resist_G = ["RR", "Rr"];

# set a threshold after which we do not bother to calculate spread
threshold = 0.95

# set up the evaluation points for quantitative resistance
g_vals = collect(lower_g : dg : upper_g);   

# The set up the set of seeds added after the burnin period
# Four different source populations: low mean_g and low %R, low mean_g and high %R, 
# high mean_g and low %R, high mean_g and high %R

# inject scenario
inject_sd_g = 1.0;

inj_TSR = 0.1;
num_inject = 10.0;
inject_locs = [convert(Int64, x_dim / 2)];

# set up some parameters to run the model over 
targ_sur_rr = [0.01, 0.5, 0.97]; 
off_var = [0.5, 1.0, 1.5];
g_prot = [1.0, 1.2, 1.5]; 

# create the parameter list
par_list = [];
for ov in off_var
  for ig in targ_sur_rr
    for gp in g_prot

      push!(par_list, [g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, inj_TSR,
	sur_2_g(ig, herb_effect, base_sur, gp), inject_sd_g, inject_locs, int_num_rr, 
	int_mean_g, int_sd_g, seed_sur, germ_prob, resist_G, fec_max, dd_fec, fec0, fec_cost, 
	base_sur, herb_effect, gp, pro_exposed, seed_pro_short, seed_mean_dist_short, 
	pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long, scale_pollen, 
	shape_pollen, sqrt(ov), threshold])
      
    end
  end
end

out = pmap(space_time_wrapper, par_list, batch_size = 9);
 
mat_out = vcat(out ...);
df_TSR = DataFrame(mat_out);
###########Make the labels properly, add in scen, metric, ts and location x1:x_dim to labels
names!(df_TSR, convert(Array{Symbol}, vcat(["inj_TSR", "inj_g", "seed_sur", "germ_prob", "fec_max", "dd_fec", "fec0", "fec_cost",
  "s0", "herb_effect", "g_pro", "pro_expo", "P_s_seed", "mds_seed", "P_mds_seed", "mdl_seed", "P_mdl_seed", 
  "scale_pollen", "shape_pollen", "off_sd", "scen", "metric", "ts"],  [string("x", i) for i = 1:(x_dim)])));
  
output_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output" 
cd(output_loc);
writetable("rho_Va_TSR_space.csv", df_TSR);
  
  
