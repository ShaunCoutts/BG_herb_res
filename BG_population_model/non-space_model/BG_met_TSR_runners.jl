# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

using Distributions

# one time step that up dates the population
function pop_update!(RR_landscape::Array{Float64, 2}, Rr_landscape::Array{Float64, 2}, rr_landscape::Array{Float64, 2},
  RR_ab_pop::Array{Float64, 1}, Rr_ab_pop::Array{Float64, 1}, rr_ab_pop::Array{Float64, 1}, 
  eff_pop_holder::Array{Float64, 1}, RR_eff_pop::Array{Float64, 1}, Rr_eff_pop::Array{Float64, 1}, rr_eff_pop::Array{Float64, 1},
  pollen_RR::Array{Float64, 1}, pollen_Rr::Array{Float64, 1}, pollen_rr::Array{Float64, 1},
  RR_newseed::Array{Float64, 1}, Rr_newseed::Array{Float64, 1}, rr_newseed::Array{Float64, 1}, 
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1}, 
  sur_tup::Tuple{Float64, Array{Float64, 1}}, seed_sur::Float64, g_vals::Array{Float64, 1}, resist_G::Array{String, 1}, 
  germ_prob::Float64, fec_max::Float64, dd_fec::Float64, dg::Float64, 
  herb_application::Int64, t::Int64)

  #move seeds to the next timestep, killing as we do so
  seedbank_update!(RR_landscape, RR_landscape, t, seed_sur)
  seedbank_update!(Rr_landscape, Rr_landscape, t, seed_sur)
  seedbank_update!(rr_landscape, rr_landscape, t, seed_sur)
  
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
  
  if num_at_t > 0.0 # if there are some individuals alive above ground create seeds, other wise don't bother
    
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
    RR_landscape[:, t] = RR_landscape[:, t] + RR_newseed 
    Rr_landscape[:, t] = Rr_landscape[:, t] + Rr_newseed
    rr_landscape[:, t] = rr_landscape[:, t] + rr_newseed
  
  end

  return nothing
    
end


# the first three positional arguments (int_pop_RR, int_pop_Rr, int_pop_rr) are 
# arrays of a populaiton at each location in the next three positional args, which
# are indexes of these populations in the landscape
function multi_iter(int_pop_RR::Array{Float64, 1}, int_pop_Rr::Array{Float64, 1}, int_pop_rr::Array{Float64, 1},
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1}, 
  sur_tup::Tuple{Float64, Array{Float64, 1}}, seed_sur::Float64, g_vals::Array{Float64, 1}, resist_G::Array{String, 1}, 
  num_iter::Int64, germ_prob::Float64, fec_max::Float64, dd_fec::Float64, dg::Float64, 
  herb_application::Int64)

  #set aside a chunck of memory for the landscapes for each genotype 
  RR_landscape = zeros(length(g_vals), num_iter)
  Rr_landscape = zeros(length(g_vals), num_iter)
  rr_landscape = zeros(length(g_vals), num_iter)
  
  #sets the intial population for each G at the locations specified in int_pop_x_G. 
  RR_landscape[:, 1] = deepcopy(int_pop_RR)
  Rr_landscape[:, 1] = deepcopy(int_pop_Rr)
  rr_landscape[:, 1] = deepcopy(int_pop_rr)
  
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
  
    pop_update!(RR_landscape, Rr_landscape, rr_landscape, RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      eff_pop_holder, RR_eff_pop, Rr_eff_pop, rr_eff_pop, pollen_RR, pollen_Rr, pollen_rr,
      RR_newseed, Rr_newseed, rr_newseed, g_mixing_kernel, g_mixing_index, g_effect_fec, 
      sur_tup, seed_sur, g_vals, resist_G, germ_prob, fec_max, dd_fec, dg, herb_application, t)

  end
  
  return (RR_landscape, Rr_landscape, rr_landscape)
  
end

# the first three positional arguments (int_pop_RR, int_pop_Rr, int_pop_rr) are 
# arrays of a populaiton at each location in the next three positional args, are 
# the arrays of populations over g for the seeds that are injected
function multi_iter_HSI(int_pop_RR::Array{Float64, 1}, int_pop_Rr::Array{Float64, 1}, int_pop_rr::Array{Float64, 1},
  HSI_RR::Array{Float64, 1}, HSI_Rr::Array{Float64, 1}, HSI_rr::Array{Float64, 1}, num_est::Int64, num_iter::Int64,
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1}, 
  sur_tup::Tuple{Float64, Array{Float64, 1}}, seed_sur::Float64, g_vals::Array{Float64, 1}, resist_G::Array{String, 1}, 
  germ_prob::Float64, fec_max::Float64, dd_fec::Float64, dg::Float64, 
  herb_app1::Int64, herb_app2::Int64)

  #set aside a chunck of memory for the landscapes for each genotype 
  RR_landscape = zeros(length(g_vals), num_iter + num_est)
  Rr_landscape = zeros(length(g_vals), num_iter + num_est)
  rr_landscape = zeros(length(g_vals), num_iter + num_est)
  
  #sets the intial population for each G at the locations specified in int_pop_x_G. 
  RR_landscape[:, 1] = deepcopy(int_pop_RR)
  Rr_landscape[:, 1] = deepcopy(int_pop_Rr)
  rr_landscape[:, 1] = deepcopy(int_pop_rr)
  
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
  
  # iterate through the timesteps before seed injection
  for t in 2:num_est
  
    pop_update!(RR_landscape, Rr_landscape, rr_landscape, RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      eff_pop_holder, RR_eff_pop, Rr_eff_pop, rr_eff_pop, pollen_RR, pollen_Rr, pollen_rr,
      RR_newseed, Rr_newseed, rr_newseed, g_mixing_kernel, g_mixing_index, g_effect_fec, 
      sur_tup, seed_sur, g_vals, resist_G, germ_prob, fec_max, dd_fec, dg, herb_app1, t)

  end
  
  # inject the seeds into the population
  RR_landscape[:, num_est] = RR_landscape[:, num_est] + HSI_RR
  Rr_landscape[:, num_est] = Rr_landscape[:, num_est] + HSI_Rr
  rr_landscape[:, num_est] = rr_landscape[:, num_est] + HSI_rr
  
  for t in (num_est + 1):(num_iter + num_est)
  
    pop_update!(RR_landscape, Rr_landscape, rr_landscape, RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      eff_pop_holder, RR_eff_pop, Rr_eff_pop, rr_eff_pop, pollen_RR, pollen_Rr, pollen_rr,
      RR_newseed, Rr_newseed, rr_newseed, g_mixing_kernel, g_mixing_index, g_effect_fec, 
      sur_tup, seed_sur, g_vals, resist_G, germ_prob, fec_max, dd_fec, dg, herb_app2, t)

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

function model_run(param::Array{Float64, 1}, int_g::Float64, int_sd::Float64, num_iter::Int64,  
  g_vals::Array{Float64, 1}, dg::Float64, resist_G::Array{String, 1}, herb_app::Int64)
 
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
  base_sur = param[13]
  offspring_sd = param[14]

  #intial populations at each intial location
  int_pop_RR = pdf(Normal(int_g, int_sd), g_vals) * int_num_RR
  int_pop_Rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_Rr
  int_pop_rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_rr
  
  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  
  # give the effect of herb as a function of g, make it symetrical stabilising function, centered on 0
  g_effect_fec = 1 ./ (1 + exp(-(fec0 - abs(g_vals) * fec_cost)))
  
  #set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)
  
  herb_run = multi_iter(int_pop_RR, int_pop_Rr, int_pop_rr, g_mixing_kernel, g_mixing_index, g_effect_fec, 
    sur_tup, seed_sur, g_vals, resist_G, num_iter, germ_prob, fec_max, dd_fec, dg, herb_app)
  
  return herb_run 
  
end


## TODO: IMPLEMENT A SEED INJECTION VERSION OF MULTI-ITER FOR THE TRANSLOCATION EXPERIMENTS

# version of the function that takes indidual parameter values rahter than a list, this version offers more flexibility 
# in terms of how parmeters are varied independently
function model_run_hot_seed_injection(int_num_RR::Float64, int_num_Rr::Float64, int_num_rr::Float64, 
  inj_num_RR::Float64, inj_num_Rr::Float64, inj_num_rr::Float64, int_g::Float64, int_sd::Float64, 
  inj_g::Float64, inj_sd::Float64, num_est::Int64, num_iter::Int64, herb_app1::Int64, herb_app2::Int64,
  germ_prob::Float64, fec0::Float64, fec_cost::Float64, fec_max::Float64, dd_fec::Float64, herb_effect::Float64, 
  g_prot::Float64, seed_sur::Float64, pro_exposed::Float64, base_sur::Float64, offspring_sd::Float64,
  g_vals::Array{Float64, 1}, dg::Float64, resist_G::Array{String, 1})
 
  #intial populations 
  int_pop_RR = pdf(Normal(int_g, int_sd), g_vals) * int_num_RR
  int_pop_Rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_Rr
  int_pop_rr = pdf(Normal(int_g, int_sd), g_vals) * int_num_rr
  
  # injected seeds
  HSI_RR = pdf(Normal(inj_g, inj_sd), g_vals) * inj_num_RR
  HSI_Rr = pdf(Normal(inj_g, inj_sd), g_vals) * inj_num_Rr
  HSI_rr = pdf(Normal(inj_g, inj_sd), g_vals) * inj_num_rr
  
  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  
  # give the effect of herb as a function of g, make it symetrical stabilising function, centered on 0
  g_effect_fec = 1 ./ (1 + exp(-(fec0 - abs(g_vals) * fec_cost)))
  
  #set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)
  
  herb_run = multi_iter_HSI(int_pop_RR, int_pop_Rr, int_pop_rr, HSI_RR, HSI_Rr, HSI_rr, 
    num_est, num_iter, g_mixing_kernel, g_mixing_index, g_effect_fec, sur_tup, seed_sur, 
    g_vals, resist_G, germ_prob, fec_max, dd_fec, dg, herb_app1, herb_app2)
  
  return herb_run 
  
end

# wrapper function to take the results of the herbicide run and make some populaton summaries for plotting 
# output a list of parameters, and then a set of measures over the number of time steps
function run_wrapper(param::Array{Float64, 1}, int_g::Float64, int_sd::Float64, num_iter::Int64,  
  g_vals::Array{Float64, 1}, dg::Float64, resist_G::Array{String, 1}, herb_app::Int64, 
  scen_lab::String)
  
  pop_run = model_run(param, int_g, int_sd, num_iter, g_vals, dg, resist_G, herb_app)
    
  out = Array{Any, 2}(7, length(param) + 4 + num_iter)

  # fill in parameter values and scenario
  for i in 1:size(out)[1]
  
    out[i, 1:2] = [int_g, int_sd]
    out[i, 3:(length(param) + 1)] = param[1:14]
    out[i, length(param) + 2] = scen_lab
    out[i, length(param) + 3] = param[15]
  
  end
  
  sur_pre = survival_pre_calc(param[13], g_vals, param[9], param[10], param[12])
  
  # fill in the different measures
  out[1, length(param) + 4] = "sur_rr"
  out[1, (length(param) + 5):end] = get_sur_rr(pop_run[3], param[13], param[9], param[10], 
  g_vals , dg)
  
  out[2, length(param) + 4] = "pro_R"
  out[2, (length(param) + 5):end] = get_pro_R(pop_run[1], pop_run[2], pop_run[3], dg) 
  
  out[3, length(param) + 4] = "pop_sur"
  out[3, (length(param) + 5):end] = get_pop_sur(pop_run[1], pop_run[2], pop_run[3], dg,
    sur_pre)
  
  out[4, length(param) + 4] = "pop_size"
  out[4, (length(param) + 5):end] = get_pop_size(pop_run[1], pop_run[2], pop_run[3], dg)
  
  out[5, length(param) + 4] = "ab_sur_pop"
  out[5, (length(param) + 5):end] = get_post_herb_pop(pop_run[1], pop_run[2], pop_run[3], dg,
    sur_pre, param[4])
    
  out[6, length(param) + 4] = "mean_g_rr"
  out[6, (length(param) + 5):end] = get_mean_g(pop_run[3], g_vals, dg) 
  
  out[7, length(param) + 4] = "var_g_rr"
  out[7, (length(param) + 5):end] = get_var_g(pop_run[3], g_vals, dg) 
  
  return out
  
end


# wrapper function to take the results of the herbicide run and make some populaton summaries for plotting 
# output a list of parameters, and then a set of measures over the number of time steps
function run_wrapper_hot_seed_injection(int_num_RR::Float64, int_num_Rr::Float64, int_num_rr::Float64, 
  inj_num_RR::Float64, inj_num_Rr::Float64, inj_num_rr::Float64, int_g::Float64, int_sd::Float64, 
  inj_g::Float64, inj_sd::Float64, num_est::Int64, num_iter::Int64, herb_app1::Int64, herb_app2::Int64,
  germ_prob::Float64, fec0::Float64, fec_cost::Float64, fec_max::Float64, dd_fec::Float64, herb_effect::Float64, 
  g_prot::Float64, seed_sur::Float64, pro_exposed::Float64, base_sur::Float64, offspring_sd::Float64,
  g_vals::Array{Float64, 1}, dg::Float64, resist_G::Array{String, 1})
  
  pop_run = model_run_hot_seed_injection(int_num_RR, int_num_Rr, int_num_rr, 
    inj_num_RR, inj_num_Rr, inj_num_rr, int_g, int_sd, inj_g, inj_sd, 
    num_est, num_iter, herb_app1, herb_app2, germ_prob, fec0, fec_cost, 
    fec_max, dd_fec, herb_effect, g_prot, seed_sur, pro_exposed, base_sur, 
    offspring_sd, g_vals, dg, resist_G)
    
  RR_pop = pop_run[1]
  Rr_pop = pop_run[2]
  rr_pop = pop_run[3]
    
  param = [int_num_RR, int_num_Rr, int_num_rr, inj_num_RR, inj_num_Rr, inj_num_rr, int_g, int_sd, 
    inj_g, inj_sd, herb_app1, herb_app2, germ_prob, fec0, fec_cost, fec_max, dd_fec, herb_effect, 
    g_prot, seed_sur, pro_exposed, base_sur, offspring_sd, num_est]
  
  out = Array{Any, 2}(7, length(param) + 1 + num_est + num_iter)
  
  # fill in parameter values and scenario
  for i in 1:size(out)[1]
  
    out[i, 1:length(param)] = param
  
  end
  
  sur_pre = survival_pre_calc(seed_sur, g_vals, herb_effect, g_prot, pro_exposed)
  
  # fill in the different measures
  out[1, length(param) + 1] = "sur_rr"
  out[1, (length(param) + 2):end] = get_sur_rr(rr_pop, base_sur, herb_effect, g_prot, 
    g_vals , dg)
  
  out[2, length(param) + 1] = "pro_R"
  out[2, (length(param) + 2):end] = get_pro_R(RR_pop, Rr_pop, rr_pop, dg) 
  
  out[3, length(param) + 1] = "pop_sur"
  out[3, (length(param) + 2):end] = get_pop_sur(RR_pop, Rr_pop, rr_pop, dg,
    sur_pre)
  
  out[4, length(param) + 1] = "pop_size"
  out[4, (length(param) + 2):end] = get_pop_size(RR_pop, Rr_pop, rr_pop, dg)
  
  out[5, length(param) + 1] = "ab_sur_pop"
  out[5, (length(param) + 2):end] = get_post_herb_pop(RR_pop, Rr_pop, rr_pop, dg,
    sur_pre, germ_prob)
    
  out[6, length(param) + 1] = "mean_g_rr"
  out[6, (length(param) + 2):end] = get_mean_g(rr_pop, g_vals, dg) 
  
  out[7, length(param) + 1] = "var_g_rr"
  out[7, (length(param) + 2):end] = get_var_g(rr_pop, g_vals, dg) 
  
  return out
  
end



# functions to make summaries of the distributions
# single time slice
function get_mean_g(n::Array{Float64, 1}, g_vals::Array{Float64, 1}, dg::Float64)
  
  p_n = n / (sum(n) * dg)
 
  return sum(p_n .* g_vals) * dg
    
end
# over all time steps
function get_mean_g(n::Array{Float64, 2}, g_vals::Array{Float64, 1}, dg::Float64)
  
  p_n = n ./ (sum(n, 1) * dg)
  
  return vec(sum(p_n .* g_vals, 1) * dg)
    
end

# single time slice
function get_var_g(n::Array{Float64, 1}, g_vals::Array{Float64, 1}, dg::Float64)

  p_n = n / (sum(n) * dg)
  
  mean_g = sum(p_n .* g_vals) * dg
  
  return sum(((g_vals - mean_g) .^ 2) .* p_n) * dg

end
# all time steps
function get_var_g(n::Array{Float64, 2}, g_vals::Array{Float64, 1}, dg::Float64)

  T = size(n)[2]
  out = zeros(T)
  
  tot_num = (sum(n, 1) * dg)
  p_n = zeros(length(g_vals))
  
  for t in 1:T
  
    if tot_num[t] > 0.0
    
      p_n[:] = n[:, t] / (sum(n[:, t]) * dg)
  
      mean_g = sum(p_n .* g_vals) * dg
  
      out[t] = sum(((g_vals - mean_g) .^ 2) .* p_n) * dg
      
    end
  
  end
  
  return out

end

# single time slice
function get_pro_R(RR_pop::Array{Float64, 1}, Rr_pop::Array{Float64, 1}, rr_pop::Array{Float64, 1}, dg::Float64)

  RR_tot = sum(RR_pop) * dg
  Rr_tot = sum(Rr_pop) * dg
  rr_tot = sum(rr_pop) * dg
  
  tot_pop = RR_tot + Rr_tot + rr_tot
  
  if tot_pop > 0.0
  
    return (2.0 * RR_tot + Rr_tot) / (2.0 * tot_pop)
  
  else
    
    return 0.0
  
  end

end
# all time steps
function get_pro_R(RR_pop::Array{Float64, 2}, Rr_pop::Array{Float64, 2}, rr_pop::Array{Float64, 2}, dg::Float64)

  RR_tot = sum(RR_pop, 1) * dg
  Rr_tot = sum(Rr_pop, 1) * dg
  rr_tot = sum(rr_pop, 1) * dg
  
  tot_pop = RR_tot + Rr_tot + rr_tot
  
  out = zeros(length(tot_pop))
 
  for t in 1:length(tot_pop)
    
    if tot_pop[t] > 0.0
    
      out[t] = (2 * RR_tot[t] + Rr_tot[t]) / (2.0 * tot_pop[t])
    
    end
    
  end
  
  return out

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
function get_pop_sur(RR_pop::Array{Float64, 1}, Rr_pop::Array{Float64, 1}, rr_pop::Array{Float64, 1}, dg::Float64,
  sur_tup::Tuple{Float64, Array{Float64, 1}})

  pre_herb_pop = get_pop_size(RR_pop, Rr_pop, rr_pop, dg)
  
  #apply herbicide
  RR_post = RR_pop * sur_tup[1]
  Rr_post = Rr_pop * sur_tup[1]
  rr_post = rr_pop .* sur_tup[2]
  
  post_herb_pop = get_pop_size(RR_post, Rr_post, rr_post, dg)
  
  if pre_herb_pop > 0.0
  
    return post_herb_pop / pre_herb_pop
  
  else 
  
    return 0.0
    
  end
    
end
# all time steps
function get_pop_sur(RR_pop::Array{Float64, 2}, Rr_pop::Array{Float64, 2}, rr_pop::Array{Float64, 2}, dg::Float64,
  sur_tup::Tuple{Float64, Array{Float64, 1}})

  pre_herb_pop = get_pop_size(RR_pop, Rr_pop, rr_pop, dg)
  
  #apply herbicide
  RR_post = RR_pop * sur_tup[1]
  Rr_post = Rr_pop * sur_tup[1]
  rr_post = rr_pop .* sur_tup[2]
  
  post_herb_pop = get_pop_size(RR_post, Rr_post, rr_post, dg)
  
  out = zeros(post_herb_pop)
  
  for t in 1:length(pre_herb_pop)
    
    if pre_herb_pop[t] > 0.0
    
      out[t] = post_herb_pop[t] / pre_herb_pop[t]
      
    end
  
  end
  
  return out

end

# single time step
function get_sur_rr(rr_pop::Array{Float64, 1}, base_sur::Float64, h_eff::Float64, g_pro::Float64, 
  g_vals::Array{Float64, 1} , dg::Float64)

  sur_vec = 1 ./ (1 + exp(-(base_sur - (h_eff - min(h_eff, g_vals * g_pro)))))
  num_sur = sum(rr_pop .* sur_vec) * dg
  
  num_pre = sum(rr_pop) * dg

  if num_pre > 0.0 
    
    return num_sur / num_pre

  else
  
    return 0.0
  
  end
    
end
# all time steps
function get_sur_rr(rr_pop::Array{Float64, 2}, base_sur::Float64, h_eff::Float64, g_pro::Float64, 
  g_vals::Array{Float64, 1} , dg::Float64)

  sur_vec = 1 ./ (1 + exp(-(base_sur - (h_eff - min(h_eff, g_vals * g_pro)))))
  num_sur = vec(sum(rr_pop .* sur_vec, 1) * dg)
  
  num_pre = vec(sum(rr_pop, 1) * dg)
  
  out = zeros(length(num_pre))
  
  for t in 1:length(num_pre)
  
    if num_pre[t] > 0.0
    
      out[t] = num_sur[t] / num_pre[t]
    
    end
    
  end
      
  return out

end

# single time step
function get_post_herb_pop(RR_pop::Array{Float64, 1}, Rr_pop::Array{Float64, 1}, rr_pop::Array{Float64, 1}, dg::Float64,
  sur_tup::Tuple{Float64, Array{Float64, 1}}, germ::Float64)

  RR_ab_pop = zeros(length(RR_pop))
  Rr_ab_pop = zeros(length(Rr_pop))
  rr_ab_pop = zeros(length(rr_pop))

  #germination   
  new_plants!(RR_ab_pop, RR_pop, germ)
  new_plants!(Rr_ab_pop, Rr_pop, germ)
  new_plants!(rr_ab_pop, rr_pop, germ)
  
  herb_applied = 2
  
  #above ground survival
  survival_at_t!(RR_ab_pop, resist_G, "RR", herb_applied, sur_tup)
  survival_at_t!(Rr_ab_pop, resist_G, "Rr", herb_applied, sur_tup)
  survival_at_t!(rr_ab_pop, resist_G, "rr", herb_applied, sur_tup)
  
  return get_pop_size(RR_ab_pop, Rr_ab_pop, rr_ab_pop, dg)    
  
end  
# all time steps
function get_post_herb_pop(RR_pop::Array{Float64, 2}, Rr_pop::Array{Float64, 2}, rr_pop::Array{Float64, 2}, dg::Float64,
  sur_tup::Tuple{Float64, Array{Float64, 1}}, germ::Float64)

  RR_ab_pop = zeros(size(RR_pop)[1])
  Rr_ab_pop = zeros(size(Rr_pop)[1])
  rr_ab_pop = zeros(size(rr_pop)[1])
  
  out = zeros(size(RR_pop)[2])
  
  for t in 1:size(RR_pop)[2]
  
    #germination   
    new_plants!(RR_ab_pop, RR_pop[:, t], germ)
    new_plants!(Rr_ab_pop, Rr_pop[:, t], germ)
    new_plants!(rr_ab_pop, rr_pop[:, t], germ)
    
    herb_applied = 2
    
    #above ground survival
    survival_at_t!(RR_ab_pop, resist_G, "RR", herb_applied, sur_tup)
    survival_at_t!(Rr_ab_pop, resist_G, "Rr", herb_applied, sur_tup)
    survival_at_t!(rr_ab_pop, resist_G, "rr", herb_applied, sur_tup)
    
    out[t] = get_pop_size(RR_ab_pop, Rr_ab_pop, rr_ab_pop, dg)    
  
  end
  
  return out
  
end  

# find the value of g that implies a given survival
function get_g_at_sur(targ_sur::Float64, base_sur::Float64, h_eff::Float64, g_pro::Float64)

  return (-log((1 / targ_sur) - 1) - base_sur + h_eff) / g_pro

end
  