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

  pop_sums = pop_summaries(pop_run[:, :, end], [0 , ls_size], dx, [lower_g, upper_g], dg, param_var, base_sur)
  
  return [param_fixed; param_var; pop_sums]

end

#function to call to do the filtering 
function param_filtering(num_par_comb::Int64, file_loc_out::AbstractString, out_name::AbstractString)
  #read in some of the required functions
  int_file_loc = pwd()
  # version for my computer
  #  @everywhere file_loc_func_p = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model" 
  #version for Zhus mac
  @everywhere file_loc_func_p = "/Users/shauncoutts/BG_pop_model/spatial_model"
  
  @everywhere cd(file_loc_func_p)
  @everywhere include("BG_met_TSR_space_pop_process.jl")
  @everywhere include("BG_met_TSR_space_dispersal_functions.jl")
  @everywhere include("BG_met_TSR_space_runners.jl")
  
  
  # fixed parameters This is an empty field scenario, could also have a full field scenario 
  # by making the intial lcoatins of rr at all locations. For the filter sanity check we 
  # make the intial population a small clump of seeds 
  int_pop_tot = 1.0
  landscape_size = 500.0
  dx = 1.0
  int_loc_RR = [249, 250, 251]
  int_loc_Rr = [249, 250, 251]
  int_loc_rr = [249, 250, 251]
  dg = 0.5
  int_g = 0.0 
  int_sd = 1.4142 
  lower_g = -10.0
  upper_g = 10.0
  offspring_sd = 1.0 
  num_iter = 50.0
  base_sur = 10.0 
  resist_G = ["RR", "Rr"] 
  herb_app_loc = collect(1:501)
  
  # upper and lower limits for each parameter that is varied 
  l_int_num_RR, u_int_num_RR = 0, 0
  l_int_num_Rr, u_int_num_Rr = 0.0001 * int_pop_tot, 0.2 * int_pop_tot 
  l_int_num_rr, u_int_num_rr = 0, 0 
  l_germ_prob, u_germ_prob = 0.45, 0.6
  l_fec0, u_fec0 = 5.0, 10.0
  l_fec_cost, u_fec_cost = 0.0, 1.0
  l_fec_max, u_fec_max = 30.0, 300.0 
  l_dd_fec, u_dd_fec = 0.001, 0.15
  l_herb_effect, u_herb_effect = 2.0, 3.0 
  l_g_prot, u_g_prot = 0.1, 2.0
  l_seed_sur, u_seed_sur = 0.22, 0.79
  l_pro_exposed, u_pro_exposed = 0.5, 1
  l_scale_pollen, u_scale_pollen = 25.6, 38.4 
  l_shape_pollen, u_shape_pollen = 2.656, 3.984 
  l_seed_pro_short, u_seed_pro_short = 0.384, 0.576
  l_seed_mean_dist_short, u_seed_mean_dist_short = 0.464, 0.696 
  l_pro_seeds_to_mean_short, u_pro_seeds_to_mean_short = 0.352, 0.528 
  l_seed_mean_dist_long, u_seed_mean_dist_long = 1.32, 1.98 
  l_pro_seeds_to_mean_long, u_pro_seeds_to_mean_long = 0.312, 0.468 
  
  #LHS of the parameter space 
  pars0 = BlackBoxOptim.Utils.latin_hypercube_sampling(
    [l_int_num_RR, l_int_num_Rr, l_int_num_rr, l_germ_prob, l_fec0, l_fec_cost, l_fec_max, 
      l_dd_fec, l_herb_effect, l_g_prot, l_seed_sur, l_pro_exposed, l_scale_pollen, l_shape_pollen, 
      l_seed_pro_short, l_seed_mean_dist_short, l_pro_seeds_to_mean_short, l_seed_mean_dist_long,
      l_pro_seeds_to_mean_long], 
    [u_int_num_RR, u_int_num_Rr, u_int_num_rr, u_germ_prob, u_fec0, u_fec_cost, u_fec_max, 
      u_dd_fec, u_herb_effect, u_g_prot, u_seed_sur, u_pro_exposed, u_scale_pollen, u_shape_pollen, 
      u_seed_pro_short, u_seed_mean_dist_short, u_pro_seeds_to_mean_short, u_seed_mean_dist_long,
      u_pro_seeds_to_mean_long],
    num_par_comb)
  #rescal a few of these as they are dependent on other parameters
  pars0[3, :] = int_pop_tot - pars0[2, :] #make the intial rr population every thing that is not Rr 
  pars0[9, :] *= base_sur # scale herb effect to s0
  pars0[10, :] = pars0[10, :] .* pars0[9, :] # scale the protective efect of g to herb_effect 
  pars0[6, :] = pars0[6, :] .* pars0[5, :] #scale demographic costs of resistance to fec0
  
  # cannot pass multipule arrays to map() (for good, but annoying reasons), so I will have to package 
  # all the parameters I want into a single list
  
  param_fixed = [int_pop_tot, landscape_size, dx, dg, int_g, int_sd, lower_g, upper_g, 
    offspring_sd, base_sur, num_iter]
  
  pars = []
  for i in 1:size(pars0)[2]
    push!(pars, (pars0[:, i], param_fixed, int_loc_RR, int_loc_Rr, int_loc_rr, resist_G, herb_app_loc))
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
  
  cd(int_file_loc)
 
  return nothing
end

srand(3214) #set random seed

#Shaun computer version
## run the population filter so I can call the file as nohup
# param_filtering(10, # number of parameter combinations
#   "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output", #file location for function source files 
#   "param_filtering_out.csv") #file location for output dataframe

#Zhu computer version   
param_filtering(24000, # number of parameter combinations, 2000 each proccess if I use 8 workers
  "/Users/shauncoutts/BG_pop_model", #file location for function source files 
  "param_filtering_out.csv") #file location for output dataframe

  
  
#   
# #simple runs and plotting for talks
# function pres_run()
#   cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model")
#   include("BG_met_TSR_space_pop_process.jl")
#   include("BG_met_TSR_space_dispersal_functions.jl")
#   include("BG_met_TSR_space_runners.jl")
#   cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/")
#   include("model_plotting_script.jl")
#   
#   # some constant parameters
#   int_pop_tot = 1.0;
#   landscape_size = 100.0;
#   dx = 1.0;
#   dg = 0.5;
#   int_g = 0.0 ;
#   int_sd = 1.4142;
#   lower_g = -10.0;
#   upper_g = 10.0;
#   offspring_sd = 1.0; 
#   num_iter = 40;
#   base_sur = 10.0; 
#   resist_G = ["RR", "Rr"]; 
#   herb_app_loc = collect(1:101);
#   # the varied parameters 
#   int_num_RR = 0.000001; 
#   int_num_Rr = 0.1;
#   int_num_rr = int_pop_tot - int_num_Rr;
#   germ_prob = 0.7; 
#   fec0 = 0.8;
#   fec_cost = 0.5; 
#   fec_max = 100.0; 
#   dd_fec = 0.004;
#   herb_effect = 2.0;  
#   g_prot = 1.0; 
#   seed_sur = 0.5; 
#   pro_exposed = 0.8;
#   scale_pollen = 32.0;
#   shape_pollen = 3.2;
#   seed_pro_short = 0.48; 
#   seed_mean_dist_short = 0.52;
#   pro_seeds_to_mean_short = 0.4; 
#   seed_mean_dist_long = 1.6;
#   pro_seeds_to_mean_long = 0.4;
#   
#   pars = [int_num_RR, int_num_Rr, int_num_rr, germ_prob, fec0, fec_cost, fec_max, 
#       dd_fec, herb_effect, g_prot, seed_sur, pro_exposed, scale_pollen, shape_pollen, 
#       seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long,
#       pro_seeds_to_mean_long];  
#   pars[3] = int_pop_tot - pars[2] - pars[1]; #make the intial rr population every thing that is not Rr 
#   pars[9] *= base_sur; # scale herb effect to s0
#   pars[10] = pars[10] * pars[9]; # scale the protective efect of g to herb_effect 
#   pars[6] = pars[6] * pars[5]; #scale demographic costs of resistance to fec0
#  
#   # first set up a run into an empty field under herbicide 
#   int_loc_RR = [49, 50, 51];
#   int_loc_Rr = [49, 50, 51];
#   int_loc_rr = [49, 50, 51];
#  
#   inv_empty_ls = model_run(pars, int_loc_RR, int_loc_Rr, int_loc_rr, int_g, int_sd, num_iter, landscape_size, dx, 
#     lower_g, upper_g, offspring_sd, dg, base_sur, resist_G, herb_app_loc);
#     
#   spatial_plot(inv_empty_ls; plot_herb = true, 
#     file_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/",
#     file_out_name = "inv_empty_ls.png")
# 
#   # invasion into a full landscape
#   int_loc_RR = [49, 50, 51];
#   int_loc_Rr = [49, 50, 51];
#   int_loc_rr = collect(1:100);
#   int_pop_tot = 10.0;
#   int_num_RR = 0.000001; 
#   int_num_Rr = 0.1;
#   int_num_rr = int_pop_tot - int_num_Rr;
#   
#   pars = [int_num_RR, int_num_Rr, int_num_rr, germ_prob, fec0, fec_cost, fec_max, 
#       dd_fec, herb_effect, g_prot, seed_sur, pro_exposed, scale_pollen, shape_pollen, 
#       seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long,
#       pro_seeds_to_mean_long];  
#   pars[3] = int_pop_tot - pars[2] - pars[1]; #make the intial rr population every thing that is not Rr 
#   pars[9] *= base_sur; # scale herb effect to s0
#   pars[10] = pars[10] * pars[9]; # scale the protective efect of g to herb_effect 
#   pars[6] = pars[6] * pars[5]; #scale demographic costs of resistance to fec0
#  
#   inv_full_ls = model_run(pars, int_loc_RR, int_loc_Rr, int_loc_rr, int_g, int_sd, num_iter, landscape_size, dx, 
#     lower_g, upper_g, offspring_sd, dg, base_sur, resist_G, herb_app_loc);
#     
#   spatial_plot(inv_full_ls; plot_herb = true, 
#     file_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/",
#     file_out_name = "inv_full_ls.png")
# 
# #plot a simple blank plot with no information plotted on it 
# dummy_data_block = (ones(size(inv_full_ls[1])[1]), zeros(size(inv_full_ls[2])), zeros(size(inv_full_ls[3]))); 
# dummy_data_block[2][:, 1, :] = 1.0;
# dummy_data_block[3][:, 1, :] = 1.0;
# 
# spatial_plot(dummy_data_block, plot_herb = true, 
#   file_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/",
#   file_out_name = "space_time_blank.png")
# 
# 
# end
# 
# 



