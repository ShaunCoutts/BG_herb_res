# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

using Distributions

# the first three positional arguments (int_pop_RR, int_pop_Rr, int_pop_rr) are 
# arrays of a populaiton at each location in the next three positional args, which
# are indexes of these populations in the landscape
function multi_iter_1D(int_pop_RR::Array{Float64, 1}, int_pop_Rr::Array{Float64, 1}, 
  int_pop_rr::Array{Float64, 1}, int_loc_RR::Array{Int64, 1}, int_loc_Rr::Array{Int64, 1},
  int_loc_rr::Array{Int64, 1}, g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2},
  g_effect_fec::Array{Float64, 1}, sur_tup::Tuple{Float64, Array{Float64, 1}}, seed_sur::Float64, 
  g_vals::Array{Float64, 1}, resist_G::Array{String, 1}, num_iter::Int64, landscape_size::Int64, 
  dx::Float64, germ_prob::Float64, fec_max::Float64, dd_fec::Float64, dg::Float64, 
  herb_application::Array{Int64, 1}, scale_pollen::Float64, shape_pollen::Float64, seed_pro_short::Float64, 
  seed_mean_dist_short::Float64, pro_seeds_to_mean_short::Float64, seed_mean_dist_long::Float64, 
  pro_seeds_to_mean_long::Float64)

   
  # pre-calc the dispersal kernels for pollen and seed 
  seed_disp_mat_1D = zeros(convert(Int32, (landscape_size / dx) + 1), 
    convert(Int32, (landscape_size / dx) + 1))
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, seed_mean_dist_short, 
    pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long)
  
  pollen_disp_mat = zeros(convert(Int32, (landscape_size / dx) + 1), 
    convert(Int32, (landscape_size / dx) + 1))
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

  ## create the RR_eff_pop and do the pre calc effect of resist costs
  RR_eff_pop = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  Rr_eff_pop = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  rr_eff_pop = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  eff_pop_holder = zeros(length(g_vals), size(seed_disp_mat_1D)[1])
  
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
    
    ## make effective pop by mult actual pop by resist and density effects
    num_at_x = (sum(RR_ab_pop, 1) + sum(Rr_ab_pop, 1) + sum(rr_ab_pop, 1)) * dg
    eff_pop_holder[:, :] = density_effect(dd_fec, num_at_x) .* g_effect_fec
    RR_eff_pop[:, :] = RR_ab_pop .* eff_pop_holder   
    Rr_eff_pop[:, :] = Rr_ab_pop .* eff_pop_holder   
    rr_eff_pop[:, :] = rr_ab_pop .* eff_pop_holder   
    
    ## pollen at arrived each location for each g from all other locations  
    pollen_RR[:, :] = RR_eff_pop * pollen_disp_mat
    pollen_Rr[:, :] = Rr_eff_pop * pollen_disp_mat
    pollen_rr[:, :] = rr_eff_pop * pollen_disp_mat
    total_pollen[:] = sum(pollen_RR, 1) + sum(pollen_Rr, 1) + sum(pollen_rr, 1)
    #normalise the pollen counts, double intergration across g, intergration across x already done is building dispersal matrix
    pollen_RR[:, :] = pollen_RR ./ transpose(total_pollen * dg)
    pollen_Rr[:, :] = pollen_Rr ./ transpose(total_pollen * dg)
    pollen_rr[:, :] = pollen_rr ./ transpose(total_pollen * dg)
    
    #create new seeds
    new_seeds_at_t_mm!(RR_newseed, Rr_newseed, rr_newseed, RR_eff_pop, Rr_eff_pop,
      rr_eff_pop, pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index,   
      fec_max, dg)
    
    # disperse the seeds
    RR_landscape[t][:, :] = RR_landscape[t]  + (RR_newseed * seed_disp_mat_1D)
    Rr_landscape[t][:, :] = Rr_landscape[t]  + (Rr_newseed * seed_disp_mat_1D)
    rr_landscape[t][:, :] = rr_landscape[t]  + (rr_newseed * seed_disp_mat_1D)
    
  end
  
  #for each one return a summary, mean g, sd, and totalnumber at each location
  #do summary here so such big arrays are not passed around
  num_loc = size(RR_landscape[1])[2]
  out_summary = zeros(9, num_loc, num_iter) #intilize the 3D array to hold the results
  
  for t in 1:num_iter
    for x in 1:num_loc
      out_summary[1:3, x, t] = dist_summary(RR_landscape[t][:, x], g_vals, dg) 
      out_summary[4:6, x, t] = dist_summary(Rr_landscape[t][:, x], g_vals, dg) 
      out_summary[7:9, x, t] = dist_summary(rr_landscape[t][:, x], g_vals, dg) 
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
  resist_G::Array{String, 1} = ["RR", "Rr"], herb_app_loc::Array{Int64, 1} = collect(1:10))
 
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
  int_pop_RR = pdf(Normal(int_g, int_sd), g_vals) * int_num_RR
  int_pop_Rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_Rr
  int_pop_rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_rr
 
  # herbicide application, one element for every location
  herb_app_0 = round(Int64, ones(convert(Int32, (landscape_size / dx) + 1)))
  herb_app = deepcopy(herb_app_0)
  herb_app[herb_app_loc] += 1 #adds one to the locations where herbicide is going to be applied 
  
  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  # give the effect of herb as a function of g, make it symetrical stabilising function, centered on 0
  g_effect_fec = exp(-(fec0 - abs(g_vals) * fec_cost))
  
  #set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)
   
  #run without herbicide
  no_herb_run = multi_iter_1D(int_pop_RR, int_pop_Rr, int_pop_rr, int_loc_RR, int_loc_Rr,
    int_loc_rr, g_mixing_kernel, g_mixing_index, g_effect_fec, sur_tup, seed_sur, g_vals, 
    resist_G, num_iter, round(Int64,landscape_size), dx, germ_prob, fec_max, dd_fec, dg, herb_app_0, 
    scale_pollen, shape_pollen, seed_pro_short, seed_mean_dist_short,pro_seeds_to_mean_short,
    seed_mean_dist_long, pro_seeds_to_mean_long)
  #run with herbicide
  herb_run = multi_iter_1D(int_pop_RR, int_pop_Rr, int_pop_rr, int_loc_RR, int_loc_Rr,
    int_loc_rr, g_mixing_kernel, g_mixing_index, g_effect_fec, sur_tup, seed_sur, g_vals, 
    resist_G, num_iter, round(Int64, landscape_size), dx, germ_prob, fec_max, dd_fec, dg, herb_app, 
    scale_pollen, shape_pollen, seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, 
    seed_mean_dist_long, pro_seeds_to_mean_long)
  
  return (param, no_herb_run, herb_run) 
end

# run the model for the filtering proccess, only need the model run under herbicide application, and there is no need to return the 
# parameters as that will be taken care of by another function
function model_run_filter(param::Array{Float64, 1}, int_loc_RR::Array{Int64, 1}, int_loc_Rr::Array{Int64, 1}, int_loc_rr::Array{Int64, 1}, 
  int_g::Float64 = 0.0, int_sd::Float64 = 1.4142, num_iter::Int64 = 10, landscape_size::Float64 = 10.0, dx::Float64 = 1.0, 
  lower_g::Float64 = -10.0, upper_g::Float64 = 10.0, offspring_sd::Float64 = 1.0, dg::Float64 = 0.5, base_sur::Float64 = 10.0, 
  resist_G::Array{String, 1} = ["RR", "Rr"], herb_app_loc::Array{Int64, 1} = collect(1:10))
 
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
  int_pop_RR = pdf(Normal(int_g, int_sd), g_vals) * int_num_RR
  int_pop_Rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_Rr
  int_pop_rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_rr
 
  # herbicide application, one element for every location
  herb_app_0 = ones(convert(Int64, (landscape_size / dx) + 1))
  herb_app = convert(Array{Int64, 1}, deepcopy(herb_app_0))
  herb_app[herb_app_loc] += 1 #adds one to the locations where herbicide is going to be applied 
  
  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
 
  # give the effect of herb resist on fecundity as a function of g
  g_effect_fec = exp(-(fec0 - abs(g_vals) * fec_cost))
  
  #set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)
   
  #run with herbicide
  herb_run = multi_iter_1D(int_pop_RR, int_pop_Rr, int_pop_rr, int_loc_RR, int_loc_Rr,
    int_loc_rr, g_mixing_kernel, g_mixing_index, g_effect_fec, sur_tup, seed_sur, g_vals, 
    resist_G, num_iter, round(Int64, landscape_size), dx, germ_prob, fec_max, dd_fec, dg, herb_app, 
    scale_pollen, shape_pollen, seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, 
    seed_mean_dist_long, pro_seeds_to_mean_long)
  
  return herb_run 
  
end

# single iteration takes a state of the population (i.e. one landscape with all TS genotypes) and progresses them
# 1 timestep foward, the first argument is the numerical next time step.
# the function takes the landscapes as references to arrays that can be modified in place (the xx_ls arrays), it all so takes a 
# set of holding arrays (xx_ab_pop, total_pollen, pollen_xx, xx_newseed, so these don't have to be built allocated every function call)

function one_step_foward!(next_t::Int64, RR_ls::Array{Float64, 3},  Rr_ls::Array{Float64, 3}, rr_ls::Array{Float64, 3}, 
  RR_ab_pop::Array{Float64, 2}, Rr_ab_pop::Array{Float64, 2}, rr_ab_pop::Array{Float64, 2}, 
  RR_eff_pop::Array{Float64, 2}, Rr_eff_pop::Array{Float64, 2}, rr_eff_pop::Array{Float64, 2}, eff_pop_holder::Array{Float64, 2},
  pollen_RR::Array{Float64, 2}, pollen_Rr::Array{Float64, 2}, pollen_rr::Array{Float64, 2}, total_pollen::Array{Float64, 1}, 
  RR_newseed::Array{Float64, 2}, Rr_newseed::Array{Float64, 2}, rr_newseed::Array{Float64, 2}, 
  dg::Float64, seed_sur::Float64, germ_prob::Float64, sur_tup::Tuple{Float64, Array{Float64, 1}}, resist_G::Array{String, 1},
   herb_application::Array{Int64, 1}, pollen_disp_mat::Array{Float64, 2}, seed_disp_mat_1D::Array{Float64, 2}, 
   fec_max::Float64, dd_fec::Float64, g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, 
   g_effect_fec::Array{Float64, 1})

  #move seeds to the next timestep, killing as we do so
  RR_ls[:, :, next_t] = RR_ls[:, :, next_t - 1] * seed_sur
  Rr_ls[:, :, next_t] = Rr_ls[:, :, next_t - 1] * seed_sur
  rr_ls[:, :, next_t] = rr_ls[:, :, next_t - 1] * seed_sur
  
  #germination   
  new_plants!(RR_ab_pop, RR_ls[:, :, next_t], germ_prob)
  new_plants!(Rr_ab_pop, Rr_ls[:, :, next_t], germ_prob)
  new_plants!(rr_ab_pop, rr_ls[:, :, next_t], germ_prob)
  
  #above ground survival
  survival_at_t!(RR_ab_pop, resist_G, "RR", herb_application, sur_tup)
  survival_at_t!(Rr_ab_pop, resist_G, "Rr", herb_application, sur_tup)
  survival_at_t!(rr_ab_pop, resist_G, "rr", herb_application, sur_tup)
  
  ## make effective pop by mult actual pop by resist and density effects
  num_at_x = (sum(RR_ab_pop, 1) + sum(Rr_ab_pop, 1) + sum(rr_ab_pop, 1)) * dg
  eff_pop_holder[:, :] = density_effect(dd_fec, num_at_x) .* g_effect_fec
  RR_eff_pop[:, :] = RR_ab_pop .* eff_pop_holder   
  Rr_eff_pop[:, :] = Rr_ab_pop .* eff_pop_holder   
  rr_eff_pop[:, :] = rr_ab_pop .* eff_pop_holder   
  
  ## pollen at arrived each location for each g from all other locations  
  pollen_RR[:, :] = RR_eff_pop * pollen_disp_mat
  pollen_Rr[:, :] = Rr_eff_pop * pollen_disp_mat
  pollen_rr[:, :] = rr_eff_pop * pollen_disp_mat
  total_pollen[:] = sum(pollen_RR, 1) + sum(pollen_Rr, 1) + sum(pollen_rr, 1)
  #normalise the pollen counts, double intergration across g, intergration across x already done in building dispersal matrix
  pollen_RR[:, :] = pollen_RR ./ transpose(total_pollen * dg)
  pollen_Rr[:, :] = pollen_Rr ./ transpose(total_pollen * dg)
  pollen_rr[:, :] = pollen_rr ./ transpose(total_pollen * dg)
  
  #create new seeds
  new_seeds_at_t_mm!(RR_newseed, Rr_newseed, rr_newseed, RR_eff_pop, Rr_eff_pop,
    rr_eff_pop, pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index,   
    fec_max, dg)
  
  # disperse the seeds
  RR_ls[:, :, next_t] = RR_ls[:, :, next_t] + (RR_newseed * seed_disp_mat_1D)
  Rr_ls[:, :, next_t] = Rr_ls[:, :, next_t] + (Rr_newseed * seed_disp_mat_1D)
  rr_ls[:, :, next_t] = rr_ls[:, :, next_t] + (rr_newseed * seed_disp_mat_1D)
  
  return nothing
  
end

# injects num_seeds into landscape_xx at locations locs, at time t_inject,
# with a given mean_g and sd_g.
function hot_seed_injection!(RR_ls_empty::Array{Float64, 3}, Rr_ls_empty::Array{Float64, 3}, 
  rr_ls_empty::Array{Float64, 3}, RR_ls_naive::Array{Float64, 3}, Rr_ls_naive::Array{Float64, 3}, 
  rr_ls_naive::Array{Float64, 3}, RR_ls_expos::Array{Float64, 3}, Rr_ls_expos::Array{Float64, 3}, 
  rr_ls_expos::Array{Float64, 3}, num_RR::Float64, num_Rr::Float64, num_rr::Float64, 
  locs::Array{Int64, 1}, t_inject::Int64, inject_mean_g::Float64, inject_sd_g::Float64, 
  g_vals::Array{Float64, 1})
  
  RR_ls_empty[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_RR
  Rr_ls_empty[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_Rr
  rr_ls_empty[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_rr
  
  RR_ls_naive[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_RR
  Rr_ls_naive[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_Rr
  rr_ls_naive[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_rr
  
  RR_ls_expos[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_RR
  Rr_ls_expos[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_Rr
  rr_ls_expos[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_rr
  
  return nothing

end


# run an entire sceanrio to produce a set of data for a plot of three measures of populaiton 
# resistance and extent over time (enough data for 9 plots, 3 measures x 3 sink pop types)
# return is a tuple of (11 x tim_steps array empty ls, 11 x tim_steps array naive ls, 11 x tim_steps array exposed ls)
function run_scene_trans(g_vals::Array{Float64, 1}, x_dim::Int64, dg::Float64, dx::Float64, 
  num_iter::Int64, burnin::Int64, num_inject::Float64, pro_R_inject::Float64, inject_mean_g::Float64, 
  inject_sd_g::Float64, inject_locs::Array{Int64, 1}, int_rr::Float64, int_mean_g::Float64,
  int_sd_g::Float64, seed_sur::Float64, germ_prob::Float64, resist_G::Array{String, 1}, 
  fec_max::Float64, dd_fec::Float64, fec0::Float64, fec_cost::Float64, base_sur::Float64, 
  herb_effect::Float64, g_prot::Float64, pro_exposed::Float64, seed_pro_short::Float64, 
  seed_mean_dist_short::Float64, pro_seeds_to_mean_short::Float64, seed_mean_dist_long::Float64, 
  pro_seeds_to_mean_long::Float64, scale_pollen::Float64, shape_pollen::Float64, offspring_sd::Float64)

  # set of temporary holdoing matricies to set aside some memory, so these don't have to be rebuilt at each iteration
  # a set of matrices to hold the above ground populations
  RR_ab_pop = zeros(length(g_vals), x_dim)
  Rr_ab_pop = zeros(length(g_vals), x_dim)
  rr_ab_pop = zeros(length(g_vals), x_dim)

  ## create the RR_eff_pop and do the pre calc effect of resist costs
  RR_eff_pop = zeros(length(g_vals), x_dim)
  Rr_eff_pop = zeros(length(g_vals), x_dim)
  rr_eff_pop = zeros(length(g_vals), x_dim)
  eff_pop_holder = zeros(length(g_vals), x_dim)
  
  # a set of matrices to hold the total amount of pollen that arrives are each location for each metabolic 
  # resitance score for each genotype
  pollen_RR = zeros(length(g_vals), x_dim)
  pollen_Rr = zeros(length(g_vals), x_dim)
  pollen_rr = zeros(length(g_vals), x_dim)
  total_pollen = zeros(x_dim)

  #set of matrices to hold the new seeds produced at each location pre dispersal 
  RR_newseed = zeros(length(g_vals), x_dim)
  Rr_newseed = zeros(length(g_vals), x_dim)
  rr_newseed = zeros(length(g_vals), x_dim)

  # set up the seed and pollen dispersal kernels 
  seed_disp_mat_1D = zeros(x_dim, x_dim)
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, seed_mean_dist_short, 
    pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long)

  pollen_disp_mat = zeros(x_dim, x_dim)
  pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx, a = scale_pollen, c = shape_pollen)

  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  
  # give the effect of herb as a function of g, make it symetrical stabilising function, centered on 0
  g_effect_fec = resist_cost_pre_calc(fec0, fec_cost, g_vals);

  # set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)

  # set up the three sink landscapes, empty, naive and exposed      
  # set aside a chunck of memory for the landscapes for each genotype or structure [g_vals, x_dim, timesteps] 
  RR_empty = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_empty = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_empty = zeros(length(g_vals), x_dim, num_iter + burnin)

  RR_naive = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_naive = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_naive = zeros(length(g_vals), x_dim, num_iter + burnin)

  RR_expos = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_expos = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_expos = zeros(length(g_vals), x_dim, num_iter + burnin)

  # intitilase the population on the full landscapes in the first time slice
  for x in 1:x_dim
    rr_naive[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_rr
    rr_expos[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_rr
  end 
  
  #herb application for navie and exposed populations
  herb_app_naive = convert(Array{Int64, 1}, ones(x_dim))
  herb_app_expos = deepcopy(herb_app_naive)
  herb_app_expos += 1 #adds one to the locations where herbicide is going to be applied 
  
  # run the burnin period on the exposed and niave pops
  for t in 2:burnin
  
    # step through the naive population, letting it develope
    one_step_foward!(t, RR_naive,  Rr_naive, rr_naive, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_naive, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
    # step through the exposed population letting it develope
    one_step_foward!(t, RR_expos,  Rr_expos, rr_expos,  
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
  end
  
  # introduce seeds from the population with different TSR and quantitative resistance
  num_RR_inject = num_inject * pro_R_inject
  num_Rr_inject = 0.0
  num_rr_inject = num_inject - num_RR_inject
  
  hot_seed_injection!(RR_empty, Rr_empty, rr_empty, RR_naive, Rr_naive, 
    rr_naive, RR_expos, Rr_expos, rr_expos, num_RR_inject, num_Rr_inject, 
    num_rr_inject, inject_locs, burnin, inject_mean_g, inject_sd_g, g_vals)
  
  # run the model for num_iter under herbicide application  
  for t in (burnin + 1):(num_iter + burnin)
  
    # step through the naive population
    one_step_foward!(t, RR_naive,  Rr_naive, rr_naive, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
    # step through the exposed population
    one_step_foward!(t, RR_expos,  Rr_expos, rr_expos, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
    # step through the empty population
    one_step_foward!(t, RR_empty,  Rr_empty, rr_empty, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
  end
  
  output_empty = pop_snapshots(RR_empty, Rr_empty, rr_empty, g_vals, dg, dx, 
    herb_effect, base_sur, g_prot)
  output_naive = pop_snapshots(RR_naive, Rr_naive, rr_naive, g_vals, dg, dx,
    herb_effect, base_sur, g_prot)
  output_expos = pop_snapshots(RR_expos, Rr_expos, rr_expos, g_vals, dg, dx,
   herb_effect, base_sur, g_prot)
  
  return (output_empty, output_naive, output_expos)
  
end

# simpler version that just injects seeds into a population once a certin quant res 
# threshold reached
function run_HSI_space(g_vals::Array{Float64, 1}, 
		       x_dim::Int64, 
		       dg::Float64, 
		       dx::Float64, 
		       num_iter::Int64, 
		       HSI_rr_thresh::Float64, 
		       HSI_Rfreq::Float64, 
		       inject_mean_g::Float64, 
		       inject_sd_g::Float64, 
		       inject_locs::Array{Int64, 1}, 
		       int_rr::Float64, 
		       int_mean_g::Float64,
		       int_sd_g::Float64, 
		       seed_sur::Float64, 
		       germ_prob::Float64, 
		       resist_G::Array{String, 1}, 
		       fec_max::Float64, 
		       dd_fec::Float64, 
		       fec0::Float64, 
		       fec_cost::Float64, 
		       base_sur::Float64, 
		       herb_effect::Float64, 
		       g_prot::Float64, 
		       pro_exposed::Float64, 
		       seed_pro_short::Float64, 
		       seed_mean_dist_short::Float64, 
		       pro_seeds_to_mean_short::Float64, 
		       seed_mean_dist_long::Float64, 
		       pro_seeds_to_mean_long::Float64, 
		       scale_pollen::Float64, 
		       shape_pollen::Float64, 
		       offspring_sd::Float64)

	# set of temporary holdoing matricies to set aside some memory, so 
	# these don't have to be rebuilt at each iteration
	# a set of matrices to hold the above ground populations
	RR_ab_pop = zeros(length(g_vals), x_dim)
	Rr_ab_pop = zeros(length(g_vals), x_dim)
	rr_ab_pop = zeros(length(g_vals), x_dim)

	## create the RR_eff_pop and do the pre calc effect of resist costs
	RR_eff_pop = zeros(length(g_vals), x_dim)
	Rr_eff_pop = zeros(length(g_vals), x_dim)
	rr_eff_pop = zeros(length(g_vals), x_dim)
	eff_pop_holder = zeros(length(g_vals), x_dim)
  
	# a set of matrices to hold the total amount of pollen that arrives are each location for each metabolic 
	# resitance score for each genotype
	pollen_RR = zeros(length(g_vals), x_dim)
	pollen_Rr = zeros(length(g_vals), x_dim)
	pollen_rr = zeros(length(g_vals), x_dim)
	total_pollen = zeros(x_dim)

	#set of matrices to hold the new seeds produced at each location pre dispersal 
	RR_newseed = zeros(length(g_vals), x_dim)
	Rr_newseed = zeros(length(g_vals), x_dim)
	rr_newseed = zeros(length(g_vals), x_dim)

	# set up the seed and pollen dispersal kernels 
	seed_disp_mat_1D = zeros(x_dim, x_dim)
	seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, 
		seed_mean_dist_short, pro_seeds_to_mean_short, seed_mean_dist_long, 
		pro_seeds_to_mean_long)

	pollen_disp_mat = zeros(x_dim, x_dim)
	pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx, a = scale_pollen, 
		c = shape_pollen)

	# build the mixing kernel for metabolic resistance score every row is a offspring score
	# every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
	# a mean of g_m*g_p and sd of offspring sd
	g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
	fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
	g_mixing_index = generate_index_pairs(g_vals)
  
	# give the effect of herb as a function of g, make it symetrical stabilising function, centered on 0
	g_effect_fec = resist_cost_pre_calc(fec0, fec_cost, g_vals);

	# set up the survival vectors 
	sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)

	# set aside a chunck of memory for the landscapes for each genotype 
	RR_expos = zeros(length(g_vals), x_dim, num_iter)
	Rr_expos = zeros(length(g_vals), x_dim, num_iter)
	rr_expos = zeros(length(g_vals), x_dim, num_iter)

	# intitilase the population on the full landscapes in the first time slice
	for x in 1:x_dim

		rr_expos[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_rr
		
	end 
  
	#herb application for navie and exposed populations
	herb_app_expos = convert(Array{Int64, 1}, ones(x_dim))  
	herb_app_expos += 1 #adds one to the locations where herbicide is going to be applied 
  
	# run the burnin period on the exposed and niave pops
	int_t = 0
	for int_t in 2:num_iter
  
		# step through the exposed population letting it develope
		one_step_foward!(int_t, RR_expos,  Rr_expos, rr_expos,  
			RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
			RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
			pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
			RR_newseed, Rr_newseed, rr_newseed, 
			dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
			pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
			g_mixing_kernel, g_mixing_index, g_effect_fec)

		# check the level of quant gen if reaches the threshold, stop inital phase 
		if get_sur_rr(rr_expos[:, inject_locs[1], int_t], g_vals, base_sur,
			herb_effect, g_prot, dg) >= HSI_rr_thresh 

			break

		end

	end

	# find the number of seeds to inject at each location and put the 
	# desired number of seeds to get HSI_Rfreq at each location
	for loc in inject_locs

		# get total pop size
		pop_tot = sum(RR_expos[:, loc, int_t] + 
			      Rr_expos[:, loc, int_t] +
			      rr_expos[:, loc, int_t]) * dg

		# find the number of RR seeds that will result in frequency of HSI_Rfreq
		# if p = populaion of rr, then 
		# inj_num_RR = (HSI_R_freq * p) / (HSI_R_freq - 1)
		# assumption: HSI_R_freq < 1 (i.e. we cannot try and replace the entire population)
		# assumption: seeds are added as RR, will need slight adjustment to add Rr
		inj_num_RR = -(HSI_Rfreq * pop_tot) / (HSI_Rfreq - 1)

		# injected seeds
		RR_expos[:, loc, int_t] = RR_expos[:, loc, int_t] + 
			pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * inj_num_RR

	end

	# test the quant res was reached in time
	if int_t < (num_iter - 2)

		# now iterate foward with TSR
		for t in (int_t + 1):num_iter
  
			# step through the exposed population
			one_step_foward!(t, RR_expos,  Rr_expos, rr_expos, 
				RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
				RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
				pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
				RR_newseed, Rr_newseed, rr_newseed, 
				dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
				pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
				g_mixing_kernel, g_mixing_index, g_effect_fec)
      
		end

	end
  
	output_expos = snapshot_space_time(RR_expos,  Rr_expos, rr_expos, 
		g_vals, dg, herb_effect, base_sur , g_prot, fec0, fec_cost)
  
	return output_expos
  
end


# run an entire sceanrio to produce a set of data for a plot of three measures of populaiton 
# resistance and extent over time (enough data for 9 plots, 3 measures x 3 sink pop types)
# return is a tuple of (11 x tim_steps array empty ls, 11 x tim_steps array naive ls, 11 x tim_steps array exposed ls)
function run_scene_trans_space_exp(g_vals::Array{Float64, 1}, x_dim::Int64, dg::Float64, dx::Float64, 
  num_iter::Int64, burnin::Int64, num_inject::Float64, pro_R_inject::Float64, inject_mean_g::Float64, 
  inject_sd_g::Float64, inject_locs::Array{Int64, 1}, int_rr::Float64, int_mean_g::Float64,
  int_sd_g::Float64, seed_sur::Float64, germ_prob::Float64, resist_G::Array{String, 1}, 
  fec_max::Float64, dd_fec::Float64, fec0::Float64, fec_cost::Float64, base_sur::Float64, 
  herb_effect::Float64, g_prot::Float64, pro_exposed::Float64, seed_pro_short::Float64, 
  seed_mean_dist_short::Float64, pro_seeds_to_mean_short::Float64, seed_mean_dist_long::Float64, 
  pro_seeds_to_mean_long::Float64, scale_pollen::Float64, shape_pollen::Float64, offspring_sd::Float64)

  # set of temporary holdoing matricies to set aside some memory, so these don't have to be rebuilt at each iteration
  # a set of matrices to hold the above ground populations
  RR_ab_pop = zeros(length(g_vals), x_dim)
  Rr_ab_pop = zeros(length(g_vals), x_dim)
  rr_ab_pop = zeros(length(g_vals), x_dim)

  ## create the RR_eff_pop and do the pre calc effect of resist costs
  RR_eff_pop = zeros(length(g_vals), x_dim)
  Rr_eff_pop = zeros(length(g_vals), x_dim)
  rr_eff_pop = zeros(length(g_vals), x_dim)
  eff_pop_holder = zeros(length(g_vals), x_dim)
  
  # a set of matrices to hold the total amount of pollen that arrives are each location for each metabolic 
  # resitance score for each genotype
  pollen_RR = zeros(length(g_vals), x_dim)
  pollen_Rr = zeros(length(g_vals), x_dim)
  pollen_rr = zeros(length(g_vals), x_dim)
  total_pollen = zeros(x_dim)

  #set of matrices to hold the new seeds produced at each location pre dispersal 
  RR_newseed = zeros(length(g_vals), x_dim)
  Rr_newseed = zeros(length(g_vals), x_dim)
  rr_newseed = zeros(length(g_vals), x_dim)

  # set up the seed and pollen dispersal kernels 
  seed_disp_mat_1D = zeros(x_dim, x_dim)
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, seed_mean_dist_short, 
    pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long)

  pollen_disp_mat = zeros(x_dim, x_dim)
  pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx, a = scale_pollen, c = shape_pollen)

  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  
  # give the effect of herb as a function of g, make it symetrical stabilising function, centered on 0
  g_effect_fec = resist_cost_pre_calc(fec0, fec_cost, g_vals);

  # set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)

  # set up the three sink landscapes, empty, naive and exposed      
  # set aside a chunck of memory for the landscapes for each genotype or structure [g_vals, x_dim, timesteps] 
  RR_empty = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_empty = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_empty = zeros(length(g_vals), x_dim, num_iter + burnin)

  RR_naive = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_naive = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_naive = zeros(length(g_vals), x_dim, num_iter + burnin)

  RR_expos = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_expos = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_expos = zeros(length(g_vals), x_dim, num_iter + burnin)

  # intitilase the population on the full landscapes in the first time slice
  for x in 1:x_dim
    rr_naive[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_rr
    rr_expos[:, x, 1] = pdf(Normal(sur_2_g(0.5, herb_effect, base_sur, g_prot), int_sd_g), g_vals) * int_rr # start the exposed population with g such that survival is 0.5.
  end 
  
  #herb application for navie and exposed populations
  herb_app_naive = convert(Array{Int64, 1}, ones(x_dim))
  herb_app_expos = deepcopy(herb_app_naive)
  herb_app_expos += 1 #adds one to the locations where herbicide is going to be applied 
  
  # run the burnin period on the exposed and niave pops
  for t in 2:burnin
  
    # step through the naive population, letting it develope
    one_step_foward!(t, RR_naive,  Rr_naive, rr_naive, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_naive, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
    # step through the exposed population letting it develope
    one_step_foward!(t, RR_expos,  Rr_expos, rr_expos,  
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
  end
  
  # introduce seeds from the population with different TSR and quantitative resistance
  num_RR_inject = num_inject * pro_R_inject
  num_Rr_inject = 0.0
  num_rr_inject = num_inject - num_RR_inject
  
  hot_seed_injection!(RR_empty, Rr_empty, rr_empty, RR_naive, Rr_naive, 
    rr_naive, RR_expos, Rr_expos, rr_expos, num_RR_inject, num_Rr_inject, 
    num_rr_inject, inject_locs, burnin, inject_mean_g, inject_sd_g, g_vals)
  
  # run the model for num_iter under herbicide application  
  for t in (burnin + 1):(num_iter + burnin)
  
    # step through the naive population
    one_step_foward!(t, RR_naive,  Rr_naive, rr_naive, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
    # step through the exposed population
    one_step_foward!(t, RR_expos,  Rr_expos, rr_expos, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
    # step through the empty population
    one_step_foward!(t, RR_empty,  Rr_empty, rr_empty, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app_expos, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
  end
  
  output_empty = snapshot_space_time(RR_empty, Rr_empty, rr_empty, g_vals, dg, herb_effect, base_sur, g_prot, fec0, fec_cost)
  output_naive = snapshot_space_time(RR_naive, Rr_naive, rr_naive, g_vals, dg, herb_effect, base_sur, g_prot, fec0, fec_cost)
  output_expos = snapshot_space_time(RR_expos, Rr_expos, rr_expos, g_vals, dg, herb_effect, base_sur, g_prot, fec0, fec_cost)
    
  return (output_empty, output_naive, output_expos)
  
end

