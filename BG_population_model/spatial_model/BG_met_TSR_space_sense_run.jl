# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.
import Distributions
@everywhere using Distributions
using BlackBoxOptim
using DataFrames

# function to combine a model run and summarise the result so pmap can call the function 
@everywhere function param_tester(all_pars::Tuple) 

  # unpack the tuple for each parameters sets
  param_var = all_pars[1]
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
  pop_sums = pop_summaries(pop_run[:, :, end], [0 , ls_size], dx, [lower_g, upper_g], dg, param_var, base_sur)
  
  return [param_fixed; param_var; pop_sums]

end


# get the data as a dataframes
cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/");
passed_dat = readtable("sanity_check_pars_passed.csv");

@everywhere cd("/Users/shauncoutts/BG_pop_model/spatial_model");
@everywhere include("BG_met_TSR_space_pop_process.jl");
@everywhere include("BG_met_TSR_space_dispersal_functions.jl");
@everywhere include("BG_met_TSR_space_runners.jl");

# cannot pass multipule arrays to map() (for good, but annoying reasons), so I will have to package 
# all the parameters I want into a single list

param_fixed = transpose(convert(Array{Float64}, passed_dat[1, 2:12]))

#TODO: pulls the rest of these things in pars from passed_dat

pars = []
for i in 1:size(passed_dat)[1]
  push!(pars, (pa, param_fixed, int_loc_RR, int_loc_Rr, int_loc_rr, resist_G, herb_app_loc))
end

output = pmap(param_tester, pars)
  
df_out = DataFrame(transpose(hcat(output...)))
names!(df_out, convert(Array{Symbol, 1}, ["int_pop_tot", "ls_size", "dx", "dg", "int_g", 
  "int_sd", "lower_g", "upper_g", "offspring_sd", "base_sur", "num_iter", "int_RR", "int_Rr", "int_rr", 
  "germ_prob", "fec0", "fec_cost", "fec_max", "dd_fec", "herb_effect", "g_prot", "seed_sur", 
  "pro_exposed", "scale_pollen", "shape_pollen", "seed_pro_short", "seed_mean_dist_short", 
  "pro_seeds_to_mean_short", "seed_mean_dist_long", "pro_seeds_to_mean_long", "num_sb_RR", 
  "num_sb_Rr", "num_sb_rr", "num_sb_tot", "num_ab_tot", "num_ab_ph_tot", "mean_g_RR", 
  "mean_g_Rr", "mean_g_rr", "mean_g_pop", "pro_RR_x", "pro_Rr_x", "pro_rr_x", "pro_all_x"]))

# write the data frame to .csv file 
cd(file_loc_out)
writetable(out_name, df_out)
