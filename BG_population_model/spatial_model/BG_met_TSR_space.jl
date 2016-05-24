# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

using Distributions
using BlackBoxOptim

# function to run and call the the other functions and scripts, will eventually run the whole thing
function main_calling_function(num_par_comb::Int64)
  cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model")
  include("BG_met_TSR_space_pop_process.jl")
  include("BG_met_TSR_space_dispersal_functions.jl")
  srand(321) #set random seed
  
  #fixed parameters This is an empty field scenario, could also have a full field scenario 
  # by making the intial lcoatins of rr at all locations
  int_pop_tot = 1.0
  landscape_size = 100.0
  dx = 1.0
  int_loc_RR = [49, 50, 51]
  int_loc_Rr = [49, 50, 51]
  int_loc_rr = [49, 50, 51]
  dg = 0.5
  int_g = 0.0 
  int_sd = 1.4142 
  lower_g = -10.0
  upper_g = 10.0
  offspring_sd = 1.0 
  num_iter = 10
  base_sur = 10.0 
  resist_G = ["RR", "Rr"] 
  herb_app_loc = collect(1:101)
  
  # upper and lower limits for each parameter that is varied 
  l_int_num_RR, u_int_num_RR = 0, 0
  l_int_num_Rr, u_int_num_Rr = 0.0001 * int_pop_tot, 0.2 * int_pop_tot 
  l_int_num_rr, u_int_num_rr = 0, 0 
  l_germ_prob, u_germ_prob = 0.45, 0.6
  l_fec0, u_fec0 = 5.0, 10.0
  l_fec_cost, u_fec_cost = 0.1, 2.0
  l_fec_max, u_fec_max = 30.0, 300.0 
  l_dd_fec, u_dd_fec = 0.001, 0.1
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
  
  #turns each coloumn into a seperate array, which is the way the function takes the arguments
  pars = [pars0[:, i] for i in 1:size(pars0)[2]]
    
  #run the model run function in parallel
  output = pmap(model_run, pars, int_loc_RR, int_loc_Rr, int_loc_rr, int_g, int_sd, num_iter, landscape_size, dx, 
    lower_g, upper_g, offspring_sd, dg, base_sur, resist_G, herb_app_loc)
  
#   seed_disp_mat_2D = zeros(convert(Int32, ((landscape_size / space_res) + 1) ^ 2), 
#     convert(Int32, ((landscape_size / space_res) + 1) ^ 2))
#   
# make a burning population 

#run under no herbicide

#run with herbicide

  # set up the parallel run of the function
  intput_array = SharedArray(Float64, num_pars)
  output_array = SharedArray(Array{Array{}, A}, num_pars)
  @sync @parallel for var = range
    body
  end
  # or
  # where svd is a function that takes and element of M as an argument, so M could 
  # be an array of lists holding a set of parameters, SVD will be the function 
  # that produces 3 3D array (g by x by t) -herb, noherb for each G
  # note that 100 such runs will take up nearly 1GB, so the thousands needed to do the parameter sweeps
  # will have to saved to disk every so often, appending to the previous file 
  # look at the HDF5 package or JLD package for a file structure that can be writtne and read easily eg:
  using HDF5
  using JLD
  
  jldopen("mydata.jld", "r+") do file
    write(file, "A", A)  # alternatively, say "@write file A"
  end

  c = h5open("mydata.h5", "r") do file
    read(file, "A")
  end
  
  #Also writing to fil in the parallel execution will need some though, might be best to 
  # break the param list into chunks say 100 in size, pass each one to pmap in a loop
  # after each pmap run write the result to the file under "chunck + i"
  
  # parallel test, start julia with 3 threads that pmap will use
  # julia -p 3
  #lain latin_hypercube_sampling to make parameter ranges
  M0 = round(Int64, BlackBoxOptim.Utils.latin_hypercube_sampling([1.0, 1.0, 1.0, 1.0], [10.0, 5.0, 7.0, 5], 10));
  #turns each coloumn into a seperate array, which is the way the function takes the arguments
  M = [M0[:, i] for i in 1:size(M0)[2]]
  
  #this does not work I have to share the function with process 2 and 2
  test = pmap(fun1, M, 5.0, 11.0);
  
  
  
end

#the @everywhere macro sends the function to all threads
@everywhere function fun1(m::Array{Int64, 1}, b::Float64, c::Float64)
  mat1 = ones(m[1], m[1]) * b
  mat2 = zeros(m[2], m[2]) + c
  mat3 = rand(m[3], m[3])
  
  return (m, mat1, mat2, mat3)
end 

#simple runs and plotting for talks
function pres_run()
  cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model")
  include("BG_met_TSR_space_pop_process.jl")
  include("BG_met_TSR_space_dispersal_functions.jl")
  include("BG_met_TSR_space_runners.jl")
  cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/")
  include("model_plotting_script.jl")
  
  # some constant parameters
  int_pop_tot = 1.0;
  landscape_size = 100.0;
  dx = 1.0;
  dg = 0.5;
  int_g = 0.0 ;
  int_sd = 1.4142;
  lower_g = -10.0;
  upper_g = 10.0;
  offspring_sd = 1.0; 
  num_iter = 40;
  base_sur = 10.0; 
  resist_G = ["RR", "Rr"]; 
  herb_app_loc = collect(1:101);
  # the varied parameters 
  int_num_RR = 0.000001; 
  int_num_Rr = 0.1;
  int_num_rr = int_pop_tot - int_num_Rr;
  germ_prob = 0.7; 
  fec0 = 0.8;
  fec_cost = 0.5; 
  fec_max = 100.0; 
  dd_fec = 0.004;
  herb_effect = 2.0;  
  g_prot = 1.0; 
  seed_sur = 0.5; 
  pro_exposed = 0.8;
  scale_pollen = 32.0;
  shape_pollen = 3.2;
  seed_pro_short = 0.48; 
  seed_mean_dist_short = 0.52;
  pro_seeds_to_mean_short = 0.4; 
  seed_mean_dist_long = 1.6;
  pro_seeds_to_mean_long = 0.4;
  
  pars = [int_num_RR, int_num_Rr, int_num_rr, germ_prob, fec0, fec_cost, fec_max, 
      dd_fec, herb_effect, g_prot, seed_sur, pro_exposed, scale_pollen, shape_pollen, 
      seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long,
      pro_seeds_to_mean_long];  
  pars[3] = int_pop_tot - pars[2] - pars[1]; #make the intial rr population every thing that is not Rr 
  pars[9] *= base_sur; # scale herb effect to s0
  pars[10] = pars[10] * pars[9]; # scale the protective efect of g to herb_effect 
  pars[6] = pars[6] * pars[5]; #scale demographic costs of resistance to fec0
 
  # first set up a run into an empty field under herbicide 
  int_loc_RR = [49, 50, 51];
  int_loc_Rr = [49, 50, 51];
  int_loc_rr = [49, 50, 51];
 
  inv_empty_ls = model_run(pars, int_loc_RR, int_loc_Rr, int_loc_rr, int_g, int_sd, num_iter, landscape_size, dx, 
    lower_g, upper_g, offspring_sd, dg, base_sur, resist_G, herb_app_loc);
    
  spatial_plot(inv_empty_ls; plot_herb = true, 
    file_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/",
    file_out_name = "inv_empty_ls.png")

  # invasion into a full landscape
  int_loc_RR = [49, 50, 51];
  int_loc_Rr = [49, 50, 51];
  int_loc_rr = collect(1:100);
  int_pop_tot = 10.0;
  int_num_RR = 0.000001; 
  int_num_Rr = 0.1;
  int_num_rr = int_pop_tot - int_num_Rr;
  
  pars = [int_num_RR, int_num_Rr, int_num_rr, germ_prob, fec0, fec_cost, fec_max, 
      dd_fec, herb_effect, g_prot, seed_sur, pro_exposed, scale_pollen, shape_pollen, 
      seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long,
      pro_seeds_to_mean_long];  
  pars[3] = int_pop_tot - pars[2] - pars[1]; #make the intial rr population every thing that is not Rr 
  pars[9] *= base_sur; # scale herb effect to s0
  pars[10] = pars[10] * pars[9]; # scale the protective efect of g to herb_effect 
  pars[6] = pars[6] * pars[5]; #scale demographic costs of resistance to fec0
 
  inv_full_ls = model_run(pars, int_loc_RR, int_loc_Rr, int_loc_rr, int_g, int_sd, num_iter, landscape_size, dx, 
    lower_g, upper_g, offspring_sd, dg, base_sur, resist_G, herb_app_loc);
    
  spatial_plot(inv_full_ls; plot_herb = true, 
    file_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/",
    file_out_name = "inv_full_ls.png")

#plot a simple blank plot with no information plotted on it 
dummy_data_block = (ones(size(inv_full_ls[1])[1]), zeros(size(inv_full_ls[2])), zeros(size(inv_full_ls[3]))); 
dummy_data_block[2][:, 1, :] = 1.0;
dummy_data_block[3][:, 1, :] = 1.0;

spatial_plot(dummy_data_block, plot_herb = true, 
  file_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output/",
  file_out_name = "space_time_blank.png")


end





