# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

using Distributions

# the first three positional arguments (int_pop_RR, int_pop_Rr, int_pop_rr) are 
# arrays of a populaiton at each location in the next three positional args, which
# are indexes of these populations in the landscape
function multi_iter(int_pop_RR::Array{Float64, 1}, int_pop_Rr::Array{Float64, 1}, int_pop_rr::Array{Float64, 1},
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1}, 
  sur_tup::Tuple{Float64, Array{Float64, 1}}, seed_sur::Float64, g_vals::Array{Float64, 1}, resist_G::Array{ASCIIString, 1}, 
  num_iter::Int64, germ_prob::Float64, fec_max::Float64, dd_fec::Float64, dg::Float64, 
  herb_application::Int64)

  #set aside a chunck of memory for the landscapes for each genotype 
  RR_landscape = zeros(length(g_vals), num_iter)
  Rr_landscape = zeros(length(g_vals), num_iter)
  rr_landscape = zeros(length(g_vals), num_iter)
  
  #sets the intial population for each G at the locations specified in int_pop_x_G. 
  RR_landscape[:, 1] = int_pop_RR
  Rr_landscape[:, 1] = int_pop_Rr
  rr_landscape[:, 1] = int_pop_rr
  
  # a set of matrices to hold the above ground populations
  RR_ab_pop = zeros(length(g_vals))
  Rr_ab_pop = zeros(length(g_vals))
  rr_ab_pop = zeros(length(g_vals))

  ## create the RR_eff_pop and do the pre calc effect of resist costs
  RR_eff_pop = zeros(length(g_vals))
  Rr_eff_pop = zeros(length(g_vals))
  rr_eff_pop = zeros(length(g_vals))
  eff_pop_holder = zeros(length(g_vals))
  
  # a set of matrices to hold the total amount of pollen that arrives are each location for each metabolic 
  # resitance score for each genotype
  pollen_RR = zeros(length(g_vals))
  pollen_Rr = zeros(length(g_vals))
  pollen_rr = zeros(length(g_vals))
  total_pollen = 0.0
  
  #set of matrices to hold the new seeds produced at each location pre dispersal 
  RR_newseed = zeros(length(g_vals))
  Rr_newseed = zeros(length(g_vals))
  rr_newseed = zeros(length(g_vals))
  
  # iterate through the timesteps
  for t in 2:num_iter
    
    #move seeds to the next timestep, killing as we do so
    seedbank_update!(RR_landscape[:, t], RR_landscape[:, t - 1], seed_sur)
    seedbank_update!(Rr_landscape[:, t], Rr_landscape[:, t - 1], seed_sur)
    seedbank_update!(rr_landscape[:, t], rr_landscape[:, t - 1], seed_sur)
    
    #germination   
    new_plants!(RR_ab_pop, RR_landscape[:, t], germ_prob)
    new_plants!(Rr_ab_pop, Rr_landscape[:, t], germ_prob)
    new_plants!(rr_ab_pop, rr_landscape[:, t], germ_prob)
    
    #above ground survival
    survival_at_t!(RR_ab_pop, resist_G, "RR", herb_application, sur_tup)
    survival_at_t!(Rr_ab_pop, resist_G, "Rr", herb_application, sur_tup)
    survival_at_t!(rr_ab_pop, resist_G, "rr", herb_application, sur_tup)
    
    ## make effective pop by mult actual pop by resist and density effects
    num_at_t = (sum(RR_ab_pop) + sum(Rr_ab_pop) + sum(rr_ab_pop)) * dg
    eff_pop_holder[:] = density_effect(dd_fec, num_at_t) * g_effect_fec
    RR_eff_pop[:] = RR_ab_pop .* eff_pop_holder   
    Rr_eff_pop[:] = Rr_ab_pop .* eff_pop_holder   
    rr_eff_pop[:] = rr_ab_pop .* eff_pop_holder   
    
    ## pollen for each g, normalise the pollen counts  
    total_pollen = (sum(RR_eff_pop) + sum(Rr_eff_pop) + sum(rr_eff_pop)) * dg
    pollen_RR[:] = RR_eff_pop / total_pollen
    pollen_Rr[:] = Rr_eff_pop / total_pollen
    pollen_rr[:] = rr_eff_pop / total_pollen
    
    #create new seeds
    new_seeds_at_t!(RR_newseed, Rr_newseed, rr_newseed, RR_eff_pop, Rr_eff_pop,
      rr_eff_pop, pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index,   
      fec_max, dg)
    
    # add new seeds to the seed bank
    RR_landscape[:, t] = RR_landscape[:, t]  + RR_newseed 
    Rr_landscape[:, t] = Rr_landscape[:, t]  + Rr_newseed
    rr_landscape[:, t] = rr_landscape[:, t]  + rr_newseed
    
  end
  
  return (RR_landscape, Rr_landscape, rr_landscape)
  
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
  resist_G::Array{ASCIIString, 1} = ["RR", "Rr"])
 
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
  
  # set up evaluation points on g
  g_vals = collect(lower_g : dg : upper_g)
 
  #intial populations at each intial location
  int_pop_RR = pdf(Normal(int_g, int_sd), g_vals) * int_num_RR
  int_pop_Rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_Rr
  int_pop_rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_rr
 
  # herbicide application, one element for every location
  herb_app_0 = round(Int64, ones(convert(Int32, (landscape_size / dx) + 1)))
  herb_app = deepcopy(herb_app_0)
  herb_app = 2 #adds one to the locations where herbicide is going to be applied 
  
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
  
  herb_run = multi_iter_1D(int_pop_RR, int_pop_Rr, int_pop_rr, int_loc_RR, int_loc_Rr,
    int_loc_rr, g_mixing_kernel, g_mixing_index, g_effect_fec, sur_tup, seed_sur, g_vals, 
    resist_G, num_iter, round(Int64, landscape_size), dx, germ_prob, fec_max, dd_fec, dg, herb_app, 
    scale_pollen, shape_pollen, seed_pro_short, seed_mean_dist_short, pro_seeds_to_mean_short, 
    seed_mean_dist_long, pro_seeds_to_mean_long)
  
  return herb_run 
  
end

# functions to make summaries of the distributions
# single time slice
function get_mean_g(n::Array{Float64, 1}, g_vals::Array{Float64, 1}, dg::Float64)
  
  return sum(n .* g_vals) * dg
    
end
# over all time steps
function get_mean_g(n::Array{Float64, 2}, g_vals::Array{Float64, 1}, dg::Float64)
  
  return vec(sum(n .* g_vals, 1) * dg)
    
end

# single time slice
function get_var_g(n::Array{Float64, 1}, g_vals::Array{Float64, 1}, dg::Float64)

  mean_g = get_mean_g(n, g_vals, dg)
  
  return sum(((g_vals - mean_g) .^ 2) .* n) * dg

end
# all time steps
function get_var_g(n::Array{Float64, 2}, g_vals::Array{Float64, 1}, dg::Float64)

  mean_g = get_mean_g(n, g_vals, dg)
  
  T = length(mean_g)
  out = zeros(T)
  for t in 1:T
  
    out[t] = sum(((g_vals - mean_g[t]) .^ 2) .* n[:, t]) * dg
  
  end
  
  return out

end

# single time slice
function get_pro_R(RR_pop::Array{Float64, 1}, Rr_pop::Array{Float64, 1}, rr_pop::Array{Float64, 1}, dg::Float64)

  RR_tot = sum(RR_pop) * dg
  Rr_tot = sum(Rr_pop) * dg
  rr_tot = sum(rr_pop) * dg
  
  return (2 * RR_tot + Rr_tot) / (RR_tot + Rr_tot + rr_tot)

end
# all time steps
function get_pro_R(RR_pop::Array{Float64, 2}, Rr_pop::Array{Float64, 2}, rr_pop::Array{Float64, 2}, dg::Float64)

  RR_tot = sum(RR_pop, 1) * dg
  Rr_tot = sum(Rr_pop, 1) * dg
  rr_tot = sum(rr_pop, 1) * dg
  
  return vec((2 * RR_tot + Rr_tot) ./ (RR_tot + Rr_tot + rr_tot))

end

# single time slice
function get_pop_size(RR_pop::Array{Float64, 1}, Rr_pop::Array{Float64, 1}, rr_pop::Array{Float64, 1}, dg::Float64)

  return (sum(RR_pop) + sum(Rr_pop) + sum(rr_pop)) * dg
  
end
# all time steps
function get_pop_size(RR_pop::Array{Float64, 2}, Rr_pop::Array{Float64, 2}, rr_pop::Array{Float64, 2}, dg::Float64)

  return vec((sum(RR_pop, 1) + sum(Rr_pop, 1) + sum(rr_pop, 1)) * dg)
  
end

# single time step
function get_pop_sur(RR_pop::Array{Float64, 1}, Rr_pop::Array{Float64, 1}, rr_pop::Array{Float64, 1}, dg::Float64
  sur_tup::Tuple{Float64, Array{Float64, 1}})

  pre_herb_pop = get_pop_size(RR_pop, Rr_pop, rr_pop, dg)
  
  #apply herbicide
  RR_post = RR_pop * sur_tup[1]
  Rr_post = Rr_pop * sur_tup[1]
  rr_post = rr_pop .* sur_tup[2]
  
  post_herb_pop = get_pop_size(RR_post, Rr_post, rr_post, dg)
  
  return post_herb_pop / pre_herb_pop

end
# all time steps
function get_pop_sur(RR_pop::Array{Float64, 2}, Rr_pop::Array{Float64, 2}, rr_pop::Array{Float64, 2}, dg::Float64
  sur_tup::Tuple{Float64, Array{Float64, 1}})

  pre_herb_pop = get_pop_size(RR_pop, Rr_pop, rr_pop, dg)
  
  #apply herbicide
  RR_post = RR_pop * sur_tup[1]
  Rr_post = Rr_pop * sur_tup[1]
  rr_post = rr_pop .* sur_tup[2]
  
  post_herb_pop = get_pop_size(RR_post, Rr_post, rr_post, dg)
  
  return vec(post_herb_pop ./ pre_herb_pop)

end



