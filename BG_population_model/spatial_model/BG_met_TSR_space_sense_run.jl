# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.
import Distributions
@everywhere using Distributions
using DataFrames

# function to combine a model run and summarise the result so pmap can call the function 
@everywhere function sense_tester(all_pars::Tuple) 

  # unpack the tuple for each parameters sets
  param_var = vec(all_pars[1])
  param_fixed = all_pars[2]
  int_loc_RR = all_pars[3]
  int_loc_Rr = all_pars[4]
  int_loc_rr = all_pars[5]
  resist_G = all_pars[6]
  herb_app_loc = all_pars[7]
 
  # unpar the fixed parameter array further
  ls_size = param_fixed[2] 
  dx = param_fixed[3] 
  dg = param_fixed[4] 
  int_g = param_fixed[5] 
  int_sd = param_fixed[6]
  lower_g = param_fixed[7] 
  upper_g = param_fixed[8] 
  offspring_sd = param_fixed[9] 
  base_sur = param_fixed[10]
  num_iter = convert(Int64, param_fixed[11])
 
  pop_run = model_run_filter(param_var, int_loc_RR, int_loc_Rr, int_loc_rr, int_g, int_sd, num_iter, 
    ls_size, dx, lower_g, upper_g, offspring_sd, dg, base_sur, resist_G, herb_app_loc);

    
  #TODO: measure the different things from population, write new function to measure these
  pop_sums = sense_summaries(pop_run, [0 , ls_size], dx, [lower_g, upper_g], dg, param_var, base_sur)
  
  return [param_fixed; param_var; pop_sums]

end

# takes the model run object produced by model_run_filter() to produce the metrics of population spread
# for the sensitivity anlaysis
@everywhere function sense_summaries(pop_obj::Array{Float64, 3}, ls_ext::Array{Float64, 1}, 
  dx::Float64, g_ext::Array{Float64, 1}, dg::Float64, param_var::Array{Float64, 1}, 
  base_sur::Float64)
  
  # extract some parameters that are useful
  germ_prob = param_var[4]
  herb_effect = param_var[9]
  g_prot = param_var[10]
  pro_exposed = param_var[12]
    
  #can work out R_50 out side the loop by taking slices of the array and looking for first index
  tot_R = 2 * sum(pop_obj[3, :, :], 2) * dx + sum(pop_obj[6, :, :], 2) * dx 
  tot_r = 2 * sum(pop_obj[9, :, :], 2) * dx + sum(pop_obj[6, :, :], 2) * dx 
  R_at_t = tot_R ./ (tot_R + tot_r) 
  R_50 = findfirst(R_at_t .>= 0.5)
  
  #spread indices
  num_locs = size(pop_obj)[2]
  pop_at_x = sum(pop_obj[[3, 6, 9], :, :], 1) #combine the populations of RR, Rr and rr
  occ_x_all = sum(pop_at_x .> 1, 2) ./ num_locs #proportion of locations that have at least 1 invidivual at each time step
  occ_x = occ_x_all[occ_x_all .< 0.9]
  spread_rate = ifelse(length(occ_x) == 0, NaN, mean(occ_x[2:end] - occ_x[1:(end - 1)]))
  
  # in order to calculate several of the population summaries I need to reconstruct the population
  landscape = collect(ls_ext[1] : dx : ls_ext[2])
  g_vals = collect(g_ext[1] : dg : g_ext[2])
  pop_sb_RR = zeros(length(g_vals), length(landscape)) # blank landscape seedbank
  pop_sb_Rr = zeros(length(g_vals), length(landscape)) 
  pop_sb_rr = zeros(length(g_vals), length(landscape)) 
  # set up the vectors to hold the measures for each time step as most the metrics revolove 
  # around measuring the populations over time.
  num_ts = size(pop_obj)[3]
  tot_rr_at_t = zeros(num_ts)
  tot_pop_at_t = zeros(num_ts)
  
  for t in 1:num_ts
    #creat the seed bank and above ground numbers both pre and post herbicide  
    for x in 1:length(landscape)
      pop_sb_RR[:, x] = pdf(Normal(pop_obj[1, x, t], pop_obj[2, x, t]), g_vals) * pop_obj[3, x, t]
      pop_sb_Rr[:, x] = pdf(Normal(pop_obj[4, x, t], pop_obj[5, x, t]), g_vals) * pop_obj[6, x, t]
      pop_sb_rr[:, x] = pdf(Normal(pop_obj[7, x, t], pop_obj[8, x, t]), g_vals) * pop_obj[9, x, t]
    end
    # create the above ground populations, the post herbcide matricies also creatged to be mutated later
    pop_ab_ph_RR = pop_sb_RR * germ_prob
    pop_ab_ph_Rr = pop_sb_Rr * germ_prob
    pop_ab_ph_rr = pop_sb_rr * germ_prob 
    
    #do one iteration of survival to ge the number of post herbicide above ground plants
    sur_pre_calc = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)
    survival_at_t!(pop_ab_ph_RR, ["RR", "Rr"], "RR", convert(Array{Int64, 1}, ones(length(landscape)) + 1), sur_pre_calc)
    survival_at_t!(pop_ab_ph_Rr, ["RR", "Rr"], "Rr", convert(Array{Int64, 1}, ones(length(landscape)) + 1), sur_pre_calc)
    survival_at_t!(pop_ab_ph_rr, ["RR", "Rr"], "rr", convert(Array{Int64, 1}, ones(length(landscape)) + 1), sur_pre_calc)
    
    tot_pop_at_t[t] = (sum(pop_ab_ph_RR) + sum(pop_ab_ph_Rr) + sum(pop_ab_ph_rr)) * dg * dx 
    tot_rr_at_t[t] = sum(pop_ab_ph_rr) * dg * dx
  end
  
  pro_rr_sur = sum(tot_rr_at_t) / sum(tot_pop_at_t)

  return [tot_pop_at_t[end], occ_x_all[end], R_50, spread_rate, pro_rr_sur]
  
end
  
# get the data as a dataframes
#cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/");
cd("/Users/shauncoutts/BG_pop_model/spatial_model");
passed_dat = readtable("sanity_check_pars_passed.csv");

@everywhere cd("/Users/shauncoutts/BG_pop_model/spatial_model");
#@everywhere cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model")
@everywhere include("BG_met_TSR_space_pop_process.jl");
@everywhere include("BG_met_TSR_space_dispersal_functions.jl");
@everywhere include("BG_met_TSR_space_runners.jl");

# cannot pass multipule arrays to map() (for good, but annoying reasons), so I will have to package 
# all the parameters I want into a single list

param_fixed = transpose(convert(Array{Float64}, passed_dat[1, 2:12]));
resist_G = ["RR", "Rr"]; 
herb_app_loc = collect(1:501);
int_loc_RR = [249, 250, 251]
int_loc_Rr = [249, 250, 251]
int_loc_rr = [249, 250, 251]


pars = []
for i in 1:size(passed_dat)[1]
push!(pars, (transpose(convert(Array{Float64}, passed_dat[i, 13:31])), param_fixed, 
    int_loc_RR, int_loc_Rr, int_loc_rr, resist_G, herb_app_loc))
end

output = pmap(sense_tester, pars[1:9]);
  
df_out = DataFrame(transpose(hcat(output...)));
names!(df_out, convert(Array{Symbol, 1}, ["int_pop_tot", "ls_size", "dx", "dg", "int_g", 
  "int_sd", "lower_g", "upper_g", "offspring_sd", "base_sur", "num_iter", "int_RR", "int_Rr", "int_rr", 
  "germ_prob", "fec0", "fec_cost", "fec_max", "dd_fec", "herb_effect", "g_prot", "seed_sur", 
  "pro_exposed", "scale_pollen", "shape_pollen", "seed_pro_short", "seed_mean_dist_short", 
  "pro_seeds_to_mean_short", "seed_mean_dist_long", "pro_seeds_to_mean_long", "final_ab_pop", "final_x_occ", 
  "R_50", "mean_spread", "pro_rr_sur"]));

# write the data frame to .csv file 
#cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/");
cd("/Users/shauncoutts/BG_pop_model/model_output");
writetable("senes_runs.csv", df_out);
