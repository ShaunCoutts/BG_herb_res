# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

using Distributions

# the first three positional arguments (int_pop_RR, int_pop_Rr, int_pop_rr) are tuples of
# (intial_pop_size::Float64, intial_pop_g::Float64, int_pop_x::Array{Int8, 1})
function multi_iter_1D(int_pop_RR::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]),
  int_pop_Rr::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]),
  int_pop_rr::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]); 
  int_sd = 1.41, num_iter = 10, landscape_size = 10, dx = 1, lower_g = -10, 
  upper_g = 10, germ_prob = 0.7, offspring_sd = 1.0, fec0 = 10, fec_cost = 2,
  fec_max = 100.0, dd_fec = 0.004, dg = 0.5, base_sur = 10.0, resist_G = ["RR", "Rr"], 
  herb_application = zeros(10), herb_effect = 20.0,  g_prot = 1.0, pro_exposed = 0.8)

    
  seed_disp_mat_1D = zeros(convert(Int32, (landscape_size / dx) + 1), 
    convert(Int32, (landscape_size / dx) + 1))
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, res = dx)
  
  pollen_disp_mat = zeros(convert(Int32, (landscape_size / space_res) + 1), 
    convert(Int32, (landscape_size / space_res) + 1))
  pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx)
   
  # build the metabolic resistance evaluation points 
  g_vals = (lower_g : dg : upper_g)
 
  #set aside a chunck of memory for the landscapes for each genotype 
  RR_landscape = Array{Array{Float64, 2}}(num_iter)
  Rr_landscape = Array{Array{Float64, 2}}(num_iter)
  rr_landscape = Array{Array{Float64, 2}}(num_iter)
  # build a landscape for each genotype, for each time step 
  # these are 2D arrays of size(seed_disp_mat_1D)[1] x (upper_g - lower_g) * seed_expand 
  for t in 1:num_iter
    RR_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
    Rr_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
    rr_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  end
  
  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  g_effect_fec = ifelse(g_vals .< 0, exp(-fec0), exp(-(fec0 - g_vals * fec_cost)))
  #set the intial populations
  int_RR_dist = Normal(int_pop_RR[2], int_sd)
  int_Rr_dist = Normal(int_pop_Rr[2], int_sd)
  int_rr_dist = Normal(int_pop_rr[2], int_sd)
  #sets the intial population for each G at the locations specified in int_pop_x_G. 
  for x in int_pop_RR[3]
    RR_landscape[1][:, x] = pdf(int_RR_dist, g_vals) * int_pop_RR[1]
  end
  for x in int_pop_Rr[3]
    Rr_landscape[1][:, x] = pdf(int_Rr_dist, g_vals) * int_pop_Rr[1]
  end
  for x in int_pop_rr[3]
    rr_landscape[1][:, x] = pdf(int_rr_dist, g_vals) * int_pop_rr[1]
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
  
  
  #set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, 
    g_prot, pro_exposed)
  
  
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
    
    
    ## TODO test I need the double intergration in the pollen counts I think I do
  end

  return (RR_landscape, Rr_landscape, rr_landscape)
end

# function to run the model experiment where a single parameter set is tested under both herbicide 
# and no-herbicide
function model_run(int_pop_RR::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]),
  int_pop_Rr::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]),
  int_pop_rr::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]); int_sd_guess::Float64 = 1.41, 
  burnin_time::Int32 = 20, )

  
  
end

# function to run and call the the other functions and scripts, will eventually run the whole thing
function main_calling_function()
  cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model")
  include("BG_met_TSR_space_pop_process.jl")
  include("BG_met_TSR_space_dispersal_functions.jl")
  srand(321) #set random seed
  
  
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
  # that produces 9 3D array (g by x by t) -burnin, herb, noherb for each G 
  pmap(svd, M)
  
end 