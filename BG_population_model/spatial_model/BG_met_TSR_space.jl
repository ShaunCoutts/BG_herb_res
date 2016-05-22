# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

using Distributions
using BlackBoxOptim

# the first three positional arguments (int_pop_RR, int_pop_Rr, int_pop_rr) are 
# arrays of a populaiton at each location in the next three positional args, which
# are indexes of these populations in the landscape
function multi_iter_1D(int_pop_RR::Array{Float64, 1}, int_pop_Rr::Array{Float64, 1}, 
  int_pop_rr::Array{Float64, 1}, int_loc_RR::Array{Int64, 1}, int_loc_Rr::Array{Int64, 1}
  int_loc_rr::Array{Int64, 1}, g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2},
  g_effect_fec::Array{Float64, 1}, sur_tup::Tuple{Float64, Array{Float64, 1}}, seed_sur::Float64, 
  g_vals::Array{Float64, 1} = collect(-10: 0.5 : 10), resist_G::Array{ASCIIString, 1} = ["RR", "Rr"], 
  num_iter::Int64, landscape_size::Int64, dx::Float64, germ_prob::Float64, fec_max::Float64, 
  dd_fec::Float64, dg::Float64, herb_application::Array{Int64, 1}, scale_pollen::Float64, 
  shape_pollen::Float64, seed_pro_short::Float64, seed_mean_dist_short::Float64, 
  pro_seeds_to_mean_short::Float64, seed_mean_dist_long::Float64, pro_seeds_to_mean_long::Float64)

    
  seed_disp_mat_1D = zeros(convert(Int32, (landscape_size / dx) + 1), 
    convert(Int32, (landscape_size / dx) + 1))
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, seed_mean_dist_short, 
    pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long)
  
  pollen_disp_mat = zeros(convert(Int32, (landscape_size / space_res) + 1), 
    convert(Int32, (landscape_size / space_res) + 1))
  pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx, a = scale_pollen, c = shape_pollen)
   
  #set aside a chunck of memory for the landscapes for each genotype 
  RR_landscape = Array{Array{Float64, 2}, 1}(num_iter)
  Rr_landscape = Array{Array{Float64, 2}, 1}(num_iter)
  rr_landscape = Array{Array{Float64, 2}, 1}(num_iter)
  
  for t in 1:num_iter
    RR_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
    Rr_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
    rr_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
  end
  
  
  #sets the intial population for each G at the locations specified in int_pop_x_G. 
  for x in int_loc_RR
    RR_landscape[1][:, x] = int_pop_RR
  end
  for x in int_loc_Rr
    Rr_landscape[1][:, x] = int_pop_Rr
  end
  for x in int_loc_rr
    rr_landscape[1][:, x] = int_pop_rr
  end
  
  # a set of matrices to hold the above ground populations
  RR_ab_pop = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  Rr_ab_pop = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  rr_ab_pop = zeros(length(g_vals), size(seed_disp_mat_1D)[1])

  # a set of matrices to hold the total amount of pollen that arrives are each location for each metabolic 
  # resitance score for each genotype
  pollen_RR = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  pollen_Rr = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  pollen_rr = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  total_pollen = zeros(size(seed_disp_mat_1D)[1])
  
  #set of matrices to hold the new seeds produced at each location pre dispersal 
  RR_newseed = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  Rr_newseed = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  rr_newseed = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  
  
 # iterate through the timesteps
  for t in 2:num_iter
    
    #move seeds to the next timestep, killing as we do so
    seedbank_update!(RR_landscape[t], RR_landscape[t - 1], seed_sur)
    seedbank_update!(Rr_landscape[t], Rr_landscape[t - 1], seed_sur)
    seedbank_update!(rr_landscape[t], rr_landscape[t - 1], seed_sur)
    
    #germination   
    new_plants!(RR_ab_pop, RR_landscape[t], germ_prob)
    new_plants!(Rr_ab_pop, Rr_landscape[t], germ_prob)
    new_plants!(rr_ab_pop, rr_landscape[t], germ_prob)
    
    #above ground survival
    survival_at_t!(RR_ab_pop, resist_G, "RR", herb_application, sur_tup)
    survival_at_t!(Rr_ab_pop, resist_G, "Rr", herb_application, sur_tup)
    survival_at_t!(rr_ab_pop, resist_G, "rr", herb_application, sur_tup)
    
    ## pollen at arrived each location for each g from all other locations  
    pollen_RR[:, :] = RR_ab_pop * pollen_disp_mat
    pollen_Rr[:, :] = Rr_ab_pop * pollen_disp_mat
    pollen_rr[:, :] = rr_ab_pop * pollen_disp_mat
    total_pollen[:] = sum(pollen_RR, 1) + sum(pollen_Rr, 1) + sum(pollen_rr, 1)
    #normalise the pollen counts, double intergration across g, intergration across x already done is building dispersal matrix
    pollen_RR[:, :] = pollen_RR ./ transpose(total_pollen * dg)
    pollen_Rr[:, :] = pollen_Rr ./ transpose(total_pollen * dg)
    pollen_rr[:, :] = pollen_rr ./ transpose(total_pollen * dg)
    
    #create new seeds
    new_seeds_at_t!(RR_newseed, Rr_newseed, rr_newseed, RR_ab_pop, Rr_ab_pop,
      rr_ab_pop, pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index,   
      g_effect_fec, fec_max, dd_fec, dg)
    
    # disperse the seeds
    RR_landscape[t][:, :] = RR_landscape  + (RR_newseed * seed_disp_mat_1D)
    Rr_landscape[t][:, :] = Rr_landscape  + (Rr_newseed * seed_disp_mat_1D)
    rr_landscape[t][:, :] = rr_landscape  + (rr_newseed * seed_disp_mat_1D)
    
  end
  
  #for each one return a summary, mean g, sd, and totalnumber at each location
  #do summary here so such big arrays are not passed around
  num_loc = size(RR_landscape[1])[2]
  out_summary = zeros(9, num_loc, num_iter) #intilize the 3D array to hold the results
  
  for t in 1:num_iter
    for x in 1:num_loc
      out_summary[1:3, x, t] = dist_summary(RR_landscape[t][:, x], g_vals, dg) 
      out_summary[4:5, x, t] = dist_summary(Rr_landscape[t][:, x], g_vals, dg) 
      out_summary[6:9, x, t] = dist_summary(rr_landscape[t][:, x], g_vals, dg) 
    end
  end
  
  return out_summary
end

# function to run the model experiment where a single parameter set is tested under both herbicide 
# and no-herbicide
# don't have the burnin period. Takes time and only give the sd, which then changes with selection anyway 
# also pop sd on g will be driven down to infinity because there are alwasy some parents with high g
# just a very very small amount of them, so the offspring distribution will always be biased in favour of 
# low g individuals when there is no herbicide. What we will assume is that 0 is the aritary value of g at
# the time we start the model run, and that it has been under weak selection so pop sd is 2 * offspring variance
# = 2 * 1 = sqrt(2) = 1.4142...

function model_run(param::Array{Float64, 1}, int_loc_RR::Array{Int64, 1}, int_loc_Rr::Array{Int64, 1}, int_loc_rr::Array{Int64, 1}, 
  int_g::Float64 = 0.0, int_sd::Float64 = 1.4142, num_iter::Int64 = 10, landscape_size::Float64 = 10.0, dx::Float64 = 1.0, 
  lower_g::Float64 = -10.0, upper_g::Float64 = 10.0, offspring_sd::Float64 = 1.0, dg::Float64 = 0.5, base_sur::Float64 = 10.0, 
  resist_G::Array{ASCIIString, 1} = ["RR", "Rr"], herb_app_loc::Array{Int64, 1} = collect(1:10))
 
  # unpack the parameter vector for readability
  int_num_RR = param[1]
  int_num_Rr = param[2]
  int_num_rr = param[3]
  germ_prob = param[4]
  fec0 = param[5]
  fec_cost = param[6]
  fec_max = param[7]
  dd_fec = param[8]
  herb_effect = param[9]
  g_prot = param[10]
  seed_sur = param[11]
  pro_exposed = param[12]
  scale_pollen = param[13] 
  shape_pollen = param[14] 
  seed_pro_short = param[15] 
  seed_mean_dist_short = param[16] 
  pro_seeds_to_mean_short = param[17] 
  seed_mean_dist_long = param[18] 
  pro_seeds_to_mean_long = param[19]
  
  # set up evaluation points on g
  g_vals = collect(lower_g : dg : upper_g)
 
  #intial populations at each intial location
  int_RR_pop = pdf(Normal(int_g, int_sd), g_vals) * int_num_RR
  int_Rr_pop = pdf(Normal(int_g, int_sd), g_vals) * int_num_RR
  int_rr_pop = pdf(Normal(int_g, int_sd), g_vals) * int_num_RR
 
  # herbicide application, one element for every location
  herb_app_0 = round(Int64, ones(convert(Int32, (landscape_size / space_res) + 1)))
  herb_app = deepcopy(herb_app_0)
  herb_app[herb_app_loc] += 1 #adds one to the locations where herbicide is going to be applied 
  
  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  # give the effect of herb as a function of g
  g_effect_fec = ifelse(g_vals .< 0, exp(-fec0), exp(-(fec0 - g_vals * fec_cost)))
  
  #set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)
   
  #run without herbicide
  no_herb_run = multi_iter_1D(int_pop_RR, int_pop_Rr, int_pop_rr, int_loc_RR, int_loc_Rr
    int_loc_rr, g_mixing_kernel, g_mixing_index, g_effect_fec, sur_tup, seed_sur, g_vals, 
    resist_G, num_iter, landscape_size, dx, germ_prob, fec_max, dd_fec, dg, herb_app_0, 
    scale_pollen, shape_pollen, seed_pro_short, seed_mean_dist_short,pro_seeds_to_mean_short,
    seed_mean_dist_long, pro_seeds_to_mean_long)
  #run with herbicide
  herb_run = multi_iter_1D(int_pop_RR, int_pop_Rr, int_pop_rr, int_loc_RR, int_loc_Rr
    int_loc_rr, g_mixing_kernel, g_mixing_index, g_effect_fec, sur_tup, seed_sur, g_vals, 
    resist_G, num_iter, landscape_size, dx, germ_prob, fec_max, dd_fec, dg, herb_app, 
    scale_pollen, shape_pollen, seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, 
    seed_mean_dist_long, pro_seeds_to_mean_long)
  
  return (param, no_herb_run, herb_run) 
end

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
      l_dd_fec, `l_herb_effect, l_g_prot, l_seed_sur, l_pro_exposed, l_scale_pollen, l_shape_pollen, 
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
