# Functions to run the simulation experiments 

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

# Takes the snapshots of the population at each timestep based on the three G landscapes for a single scenario 
# returns a 12 x time_steps 2D array, first 4 rows give total pops (3 G types + total), next 4 mean_g (3 G types + total),
# nest 1 gives the survival of exposed rr indviduals under mean g
# next 3 %occ (3 G types) 
function pop_snapshots(RR_ls::Array{Float64, 3}, Rr_ls::Array{Float64, 3},
  rr_ls::Array{Float64, 3}, g_vals::Array{Float64, 1}, dg::Float64, dx::Float64, 
  herb_ef::Float64, s0::Float64, g_pro::Float64)
  
  output = zeros(12, size(RR_ls)[3]) # 11 because #G (3), mean_g for each G (3), mean_g overall, survival under mean g + total_pop, #occ each G (3) (3 + 3 + 3 + 1 + 1 + 1 = 12)  
  
  # total populations for each G
  output[1, :] = sum(RR_ls, [1, 2]) * dg * dx
  output[2, :] = sum(Rr_ls, [1, 2]) * dg * dx
  output[3, :] = sum(rr_ls, [1, 2]) * dg * dx
  # total population
  output[4, :] = sum(output[1:3, :], 1)
 
  # mean_g
  output[5, :] = vec(sum((sum(RR_ls, 2) * dx) .* g_vals, 1) * dg) ./ output[1, :]
  output[6, :] = vec(sum((sum(Rr_ls, 2) * dx) .* g_vals, 1) * dg) ./ output[2, :]
  output[7, :] = vec(sum((sum(rr_ls, 2) * dx) .* g_vals, 1) * dg) ./ output[3, :]
  # population mean_g 
  output[8, :] = vec(sum(((sum(RR_ls, 2) + sum(Rr_ls, 2) + sum(rr_ls, 2)) * dx) .* g_vals, 1) * dg) ./ output[4, :]
 
  # survival of exposed rr under mean g for rr 
  for t in 1:size(RR_ls)[3]
  
    output[9, t] = g_2_sur(output[7, t], herb_ef, s0, g_pro)
 
  end
  
  # % landscape occupied 
  output[10, :] = sum(sum(RR_ls, 1) * dg .> 1.0, 2) / size(RR_ls)[2]
  output[11, :] = sum(sum(Rr_ls, 1) * dg .> 1.0, 2) / size(Rr_ls)[2]
  output[12, :] = sum(sum(rr_ls, 1) * dg .> 1.0, 2) / size(rr_ls)[2]
  
  return output
  
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


# returns the inital mean g from a given survival of exposed indivudals 
function sur_2_g(sur::Float64, herb_ef::Float64, s0::Float64, g_pro::Float64)
  
  return - ((log((1 / sur) - 1) - herb_ef + s0) / g_pro)
  
end

# returns the survival of a given g, for an exopsed, suceptable individual 
function g_2_sur(g::Float64, herb_ef::Float64, s0::Float64, g_pro::Float64)

  return 1 / (1 + exp(-(s0 - (herb_ef - min(herb_ef, g * g_pro)))))

end

# produces a list that summaries the output_tup with a few snapshot numbers of the population in the final timestep
function trans_ex_snapshot(output_tup::Tuple{Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}},
  landscape_size::Float64, threshold::Float64, burnin::Int64)
  # set up a matrix to hold the results
  # source scenarios are in rows (empty, naive, exposed) and metrics are in coloumns
  # (%R, mean_spread_RRorRr, mean_spread_RR, mean_spread_Rr, mean_spread_rr) 
  snapshot = Array{Any, 2}(3, 7)
 
  # scenario labels
  snapshot[1, 1] = "empty"
  snapshot[2, 1] = "naive"
  snapshot[3, 1] = "expos"
  
  # % R 
  snapshot[1, 2] = get_pro_R(output_tup[1][1, end], output_tup[1][2, end], output_tup[1][3, end]) #empty_%R
  snapshot[2, 2] = get_pro_R(output_tup[2][1, end], output_tup[2][2, end], output_tup[2][3, end]) #naive_%R
  snapshot[3, 2] = get_pro_R(output_tup[3][1, end], output_tup[3][2, end], output_tup[3][3, end]) #expos_%R
  
  # survival under final g 
  for i in 1:3
    snapshot[i, 3] = output_tup[i][9, end] # sur under final g
  end
  
  # mean spread rates
  # reconstruct the number of locations occupied at each time step, rahter than proportion
  empty_occ = output_tup[1][10:12, burnin:end] * landscape_size 
  naive_occ = output_tup[2][10:12, burnin:end] * landscape_size 
  expos_occ = output_tup[3][10:12, burnin:end] * landscape_size 
  
  thres = threshold * landscape_size
  
  # do the RrorRR first, take max of Rr and RR at each time
  snapshot[1, 4] = get_mean_spread(max(empty_occ[1, :], empty_occ[2, :]), thres)
  snapshot[2, 4] = get_mean_spread(max(naive_occ[1, :], naive_occ[2, :]), thres)
  snapshot[3, 4] = get_mean_spread(max(expos_occ[1, :], expos_occ[2, :]), thres)
 
  # do the mean spread for each G
  for i in 1:3
    snapshot[1, i + 4] = get_mean_spread(empty_occ[i, :], thres)
    snapshot[2, i + 4] = get_mean_spread(naive_occ[i, :], thres)
    snapshot[3, i + 4] = get_mean_spread(expos_occ[i, :], thres)
  end
  
  return snapshot
  
end

# function to produce a sinlge plot of %R over time and survival under mean g
function get_pro_R(RR_ts::Array{Float64, 1}, Rr_ts::Array{Float64, 1}, 
  rr_ts::Array{Float64, 1})
  
  return (2 * RR_ts + Rr_ts) ./ (2 * (RR_ts + Rr_ts + rr_ts))
  
end

function get_pro_R(RR::Float64, Rr::Float64, rr::Float64)
  
  return (2 * RR + Rr) / (2 * (RR + Rr + rr))
  
end

# function to produce a sinlge plot of %R over time and survival under mean g
function get_mean_spread(ts::Array{Float64, 1}, thresh::Float64)
  
  unsat_ts = ts[ts .< thresh]
  
  if(length(unsat_ts) > 2)
  
   return occ_dif = mean(unsat_ts[2:end] - unsat_ts[1:(end - 1)]) 
  
  else
  
    return(0)
  
  end
  
end


# run an single sceanrio to produce a set of data for a plot of three measures of populaiton 
# resistance and extent over time return is x_dim x g_vals x tim_steps array. 
function run_natspread(g_vals::Array{Float64, 1}, x_dim::Int64, dg::Float64, 
  num_iter::Int64, int_rr::Array{Float64, 1}, int_Rr::Array{Float64, 1}, 
  int_RR::Array{Float64, 1}, int_mean_g::Float64, int_sd_g::Float64, seed_sur::Float64, 
  germ_prob::Float64, resist_G::Array{String, 1}, fec_max::Float64, dd_fec::Float64, 
  g_effect_fec::Array{Float64, 1}, sur_tup::Tuple{Float64, Array{Float64, 1}}, 
  offspring_sd::Float64, seed_disp_mat_1D::Array{Float64, 2}, pollen_disp_mat::Array{Float64, 2}, 
  herb_app::Array{Int64, 1})

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

  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  # give the effect of herb as a function of g, make it symetrical stabilising function, centered on 0
  
  # set up the landscapes
  # set aside a chunck of memory for the landscapes for each genotype or structure [g_vals, x_dim, timesteps] 
  RR_ls = zeros(length(g_vals), x_dim, num_iter)
  Rr_ls = zeros(length(g_vals), x_dim, num_iter)
  rr_ls = zeros(length(g_vals), x_dim, num_iter)

  # intitilase the population on the full landscapes in the first time slice
  for x in 1:x_dim
    rr_ls[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_rr[x]
    Rr_ls[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_Rr[x]
    RR_ls[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_RR[x]
  end 
    
  # run the model
  for t in 2:num_iter
  
    # step through the exposed population letting it develope
    one_step_foward!(t, RR_ls,  Rr_ls, rr_ls, 
      RR_ab_pop, Rr_ab_pop, rr_ab_pop, 
      RR_eff_pop, Rr_eff_pop, rr_eff_pop, eff_pop_holder, 
      pollen_RR, pollen_Rr, pollen_rr, total_pollen, 
      RR_newseed, Rr_newseed, rr_newseed, 
      dg, seed_sur, germ_prob, sur_tup, resist_G, herb_app, 
      pollen_disp_mat, seed_disp_mat_1D, fec_max, dd_fec, 
      g_mixing_kernel, g_mixing_index, g_effect_fec)
      
  end
  
  return (RR_ls, Rr_ls, rr_ls)
  
end

# series of functions to reduce the output of a nat spread run to a 
# 2D matrix of a metric 

function proR_time_space(natspread_out_tup, dg::Float64)
  
  out_dims = size(natspread_out_tup[1])
  #get total numbers at each location and time, summing over g
  num_RR = reshape(sum(natspread_out_tup[1], 1) * dg, out_dims[2:3])
  num_Rr = reshape(sum(natspread_out_tup[2], 1) * dg, out_dims[2:3])
  num_rr = reshape(sum(natspread_out_tup[3], 1) * dg, out_dims[2:3])
  
  # calc proportion R and turn any areas with no individuals present to 0
  pro_R = (2 * num_RR + num_Rr) ./ (2 * (num_RR + num_Rr + num_rr))
  pro_R[isnan(pro_R)] = 0.0
  
  return pro_R

end

function totnum_time_space(natspread_out_tup, dg::Float64)
  
  out_dims = size(natspread_out_tup[1])
  #get total numbers at each location and time, summing over g
  num_RR = reshape(sum(natspread_out_tup[1], 1) * dg, out_dims[2:3])
  num_Rr = reshape(sum(natspread_out_tup[2], 1) * dg, out_dims[2:3])
  num_rr = reshape(sum(natspread_out_tup[3], 1) * dg, out_dims[2:3])
  
  return num_RR + num_Rr + num_rr

end

function sur_g_time_space(natspread_out_ls::Array{Float64, 3}, dg::Float64,
  g_vals::Array{Float64, 1}, herb_ef::Float64, s0::Float64, g_pro::Float64)

  dims = size(natspread_out_ls)
  out = zeros(dims[2:3])
  for t in 1:dims[3]
    for x in 1:dims[2]
      
      num_loc = sum(natspread_out_ls[:, x, t]) * dg
      if num_loc > 0.0
	out[x, t] = g_2_sur((sum(natspread_out_ls[:, x, t] .* g_vals) * dg) ./ num_loc, 
	  herb_ef, s0, g_pro)
      else
	out[x, t] = 0.0
      end
    end
  end
  
  return out
  
end


# it is very hard to reason about the realative advantage TSR individuals have at a given point in
# space and time since it is the result of both fecundity (which are influenced by density) and 
# survival functions interacting with distribution of each TSR genotype over g.

function TSR_adv_time_space(natspread_out_tup, dg::Float64, g_vals::Array{Float64, 1},
  germ_prob::Float64, fec_max::Float64, g_effect_fec::Array{Float64, 1},  
  dd_fec::Float64, sur_tup::Tuple{Float64, Array{Float64, 1}})
 
  dims = size(natspread_out_tup[1])
  out = zeros(dims[2:3])
  
  fec_corrected = fec_max / 3
  
  for t in 1:dims[3]
    for x in 1:dims[2]

      # calculate the distribution over g that emerge for each G
      ag_RR = natspread_out_tup[1][:, x, t] * germ_prob
      ag_Rr = natspread_out_tup[2][:, x, t] * germ_prob
      ag_rr = natspread_out_tup[3][:, x, t] * germ_prob
      
      # total number of each G
      num_RR = sum(ag_RR) * dg
      num_Rr = sum(ag_Rr) * dg
      num_rr = sum(ag_rr) * dg
      num_tot = num_RR + num_Rr + num_rr
      
      # seeds produced by each G
      seeds_g = fec_corrected ./ (1 + g_effect_fec .+ dd_fec * num_tot .+ dd_fec * (g_effect_fec * num_tot))
      seeds_RR = seeds_g .* ag_RR
      seeds_Rr = seeds_g .* ag_Rr
      seeds_rr = seeds_g .* ag_rr
      
      # expected number of seeds each TSR and TSS
      tot_repo_rr = sum(seeds_rr .* sur_tup[2]) * dg
      tot_repo_R = sum((seeds_Rr + seeds_RR) * sur_tup[1]) * dg
      
      # mean number of seeds per individual both TSR and TSS
      if num_rr != 0.0
	E_rr = tot_repo_rr / num_rr
      else
	E_rr = 0.0
      end
      if (num_RR + num_Rr) != 0.0
	E_R = tot_repo_R / (num_RR + num_Rr)
      else
	E_R = 0.0
      end
      
      if E_rr != 0.0
	out[x, t] = E_R / E_rr
      else
        out[x, t] = 0.0   
      end
      
    end
  end
 
  return out
 
end



# create a function for just the realtive cost in terms of realative seed production 
function TSR_rel_fec_time_space(natspread_out_tup, dg::Float64, g_vals::Array{Float64, 1},
  germ_prob::Float64, fec_max::Float64, g_effect_fec::Array{Float64, 1},  
  dd_fec::Float64)
 
  dims = size(natspread_out_tup[1])
  out = zeros(dims[2:3])
  
  fec_corrected = fec_max / 3
  
  for t in 1:dims[3]
    for x in 1:dims[2]

      # calculate the distribution over g that emerge for each G
      ag_RR = natspread_out_tup[1][:, x, t] * germ_prob
      ag_Rr = natspread_out_tup[2][:, x, t] * germ_prob
      ag_rr = natspread_out_tup[3][:, x, t] * germ_prob
      
      # total number of each G
      num_RR = sum(ag_RR) * dg
      num_Rr = sum(ag_Rr) * dg
      num_rr = sum(ag_rr) * dg
      num_tot = num_RR + num_Rr + num_rr
      
      # seeds produced by each G
      seeds_g = fec_corrected ./ (1 + g_effect_fec .+ dd_fec * num_tot .+ dd_fec * (g_effect_fec * num_tot))
      seeds_RR = seeds_g .* ag_RR
      seeds_Rr = seeds_g .* ag_Rr
      seeds_rr = seeds_g .* ag_rr
      
      # expected number of seeds each TSR and TSS
      tot_fec_rr = sum(seeds_rr) * dg
      tot_fec_R = sum((seeds_Rr + seeds_RR)) * dg
      
      # mean number of seeds per individual both TSR and TSS
      if num_rr != 0.0
	E_rr = tot_fec_rr / num_rr
      else
	E_rr = 0.0
      end
      if (num_RR + num_Rr) != 0.0
	E_R = tot_fec_R / (num_RR + num_Rr)
      else
	E_R = 0.0
      end
      
      if E_rr != 0.0
	out[x, t] = E_R / E_rr
      else
        out[x, t] = 0.0   
      end
      
    end
  end
 
  return out
 
end

# Function that takes the output from a natspread run and calculates a area under the surface to test 
# how much total R there is over time and space, an index of how important TSR was under a given parameter 
# set

function AUS_proR_surg(natspread_out_tup, dg::Float64, dx::Float64,
  g_vals::Array{Float64, 1}, herb_ef::Float64, s0::Float64, g_pro::Float64)

  proR = proR_time_space(natspread_out_tup, dg)
  sur_rr = sur_g_time_space(natspread_out_tup[3], dg, g_vals, herb_ef, s0, g_pro)
  
  return (sum(proR .* sur_rr) * dx) / prod(size(proR))

end

function AUS_proR(natspread_out_tup, dg::Float64, dx::Float64)

  proR = proR_time_space(natspread_out_tup, dg)
  
  return (sum(proR) * dx) / prod(size(proR))

end
