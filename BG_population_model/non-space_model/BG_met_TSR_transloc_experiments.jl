# simulation experiments for the TSR and NTSR joint evolution
# to test how TSR can invade into a population with NTSR
@everywhere using DistributedArrays

using DataFrames

# need to run each source scenario three times, once in a empty landscape, one in an exposed population
# and one in a herbicide exposed population (can pre calculate all these so don't need to re-caclulate each time 
@everywhere file_loc_func_p = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/non-space_model" 
output_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output" 

@everywhere cd(file_loc_func_p)
@everywhere include("BG_met_TSR_pop_process.jl")
@everywhere include("BG_met_TSR_runners.jl")
# include("spatial_model_plotting_script.jl")
# @everywhere include("BG_met_TSR_exper_funs.jl")
  

###################################################################################################
# script to run the translocation experiments
# parameter values
upper_g = 20.0;
lower_g = -20.0;
dg = 0.5;
int_mean_g = 0.0;
int_sd_g = 1.4142;

int_num_RR = 0.0;
int_num_Rr = 0.0;
int_num_rr = 10.0; # number of intial seeds at each location for each genoptype, assume only TS susceptible
burnin = 20;
num_iter = 100;

offspring_sd = 1.0;
fec_max = 60.0;
fec0 = 4.0;
fec_cost = 0.45;
dd_fec = 0.001;
base_sur = 10.0; 
herb_effect = 16.0; 
g_prot = 1.5; 
pro_exposed = 0.8;
seed_sur = 0.45;
germ_prob = 0.52;
resist_G = ["RR", "Rr"];
herb_app = 2;

# set up the evaluation points for quantitative resistance
g_vals = collect(lower_g : dg : upper_g);   

# set a up a limited parameter to find a nice region of parameter space to work in
g_pro = [1.0, 1.5, 2.0];
Va = [0.5, 1.0, 1.5]; # addative variance, take sqrt() to standard deviation, which is the julia parameterisation

par_list = [];
rep_count = 1.0;
for gp in g_pro
  for osv in Va

    push!(par_list, [int_num_RR, int_num_Rr, int_num_rr, germ_prob, fec0, fec_cost, fec_max, 
      dd_fec, herb_effect, gp, seed_sur, pro_exposed, base_sur, sqrt(osv), rep_count]);
    
    rep_count += 1;
    
  end
end

# copy that list to all the  workers 
@eval @everywhere par_list = $par_list

# define the other needed inputs on all workers
@everywhere upper_g = 20.0;
@everywhere lower_g = -20.0;
@everywhere dg = 0.5;
@everywhere int_mean_g = 0.0;
@everywhere int_sd_g = 1.4142;
@everywhere num_iter = 100;
@everywhere resist_G = ["RR", "Rr"];
@everywhere herb_app = 2;
@everywhere g_vals = collect(lower_g : dg : upper_g);   


@time out = @DArray [run_wrapper(par_list[x], int_mean_g, int_sd_g, num_iter, g_vals, dg, 
  resist_G, herb_app, "no TSR") for x = 1:length(par_list)];


var_names = ["int_mean_g", "int_sd_g", "intRR", "intRr", "intrr", "germ_prob", "fec0", "fec_cost", 
  "fec_max", "dd_fec", "herb_effect", "g_pro", "seed_sur", "pro_exposed", "s0", "off_sd", "scen", "rep_ID", "measure"];
all_names = vcat(var_names, [string("t", i) for i = 1:num_iter]);

res_df = DataFrame(vcat(out...));
names!(res_df, convert(Array{Symbol}, all_names));

# write the parameter sweep to a .csv file so an R plotting script can be used
cd(output_loc);
writetable("rho_Va_par_sweep_nonspace.csv", res_df);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

par_vals = [int_num_RR, int_num_Rr, int_num_rr, germ_prob, fec0, fec_cost, fec_max, 
  dd_fec, herb_effect, g_prot, seed_sur, pro_exposed, base_sur, offspring_sd, 1.0];
  
  
test_out = run_wrapper(par_vals, int_mean_g, int_sd_g, num_iter, lower_g, upper_g, dg, 
  resist_G, herb_app, "no TSR")
  
  

par_list_g_pro = []
for i in 1:size(inj_scens)[1]
  push!(par_list_g_pro, [g_vals, x_dim, dg, dx, num_iter, burnin, num_inject, inj_scens[i, 2],
    inj_scens[i, 3], inject_sd_g, inject_locs, int_num_rr, int_mean_g, int_sd_g, seed_sur, germ_prob, 
    resist_G, fec_max, dd_fec, fec0, fec_cost, base_sur, herb_effect, inj_scens[i, 1], pro_exposed, 
    seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
    pro_seeds_to_mean_long, scale_pollen, shape_pollen, offspring_sd, threshold])
end

out = pmap(runner_wrapper, par_list_g_pro, batch_size = 55)



# survial curve and fec cure plots to visulise the tradeoff
function sur_fun(s0::Float64, h::Float64, gp::Float64, g_vals::Array{Float64, 1})

  return 1 ./ (1 + exp(-(s0 - (h - min(h, gp * g_vals)))))

end

plt = plot(g_vals, sur_fun(base_sur, herb_effect, g_prot, g_vals))
plt = plot!(plt, g_vals, resist_cost_pre_calc(3.0, fec_cost, g_vals))

sur = sur_fun(base_sur, herb_effect, 1.5, g_vals)
fec = resist_cost_pre_calc(fec0, fec_cost, g_vals)

sur[g_vals .== 5]
fec[g_vals .== 5]




plot(get_mean_g(test[3], g_vals, dg))
plot(get_var_g(test[3], g_vals, dg))


sur_test = survival_pre_calc(base_sur, g_vals, herb_effect, 
  g_prot, pro_exposed)
  
plot(get_post_herb_pop(test[1], test[2], test[3], dg, sur_test, germ_prob))

plot(get_pop_size(test[1], test[2], test[3], dg))



# set up a data frame to hold the results 
var_names = ["int_mean_g", "int_sd_g", "intRR", "intRr", "intrr", "germ_prob", "fec0", "fec_cost", 
  "fec_max", "dd_fec", "herb_effect", "g_pro", "seed_sur", "pro_exposed", "s0", "off_sd", "scen", "rep_ID", "measure"]
all_names = vcat(var_names, [string("t", i) for i = 1:num_iter])

first_row = vec([zeros(length(par_vals) + 2); "none"; 0; "none"; zeros(num_iter)])

res_df = DataFrame(reshape(first_row, 1, length(first_row)))
names!(res_df, convert(Array{Symbol}, all_names)) 

new_df = DataFrame(test_out)
names!(new_df, convert(Array{Symbol}, all_names))


res_df = vcat(res_df, new_df)














































































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
  
 
  
  
  
  