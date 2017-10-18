# Functions to run the simulation experiments 
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

# function that takes the raw population matricies (over g and space and time) and produces summaries of the populations at each time and space
function snapshot_space_time(RR_ls::Array{Float64, 3}, Rr_ls::Array{Float64, 3}, rr_ls::Array{Float64, 3}, 
  g_vals::Array{Float64, 1}, dg::Float64, herb_ef::Float64, s0::Float64, g_pro::Float64, fec0::Float64, 
  fec_cost::Float64)
  
  T = size(RR_ls)[3] # time steps
  ls_size = size(RR_ls)[2]
    
  m_names = ["num_RR", "num_Rr", "num_rr", "pro_R", "sur_rr", "TSR_adv"]
  
  output = Array{Any, 2}(T * length(m_names), ls_size + 2)
  
  # name the rows to show what they do and flatten out the 3D array to 2D to make it easier to handel and plot
  output[:, 1] = repeat(m_names, inner = T)
  output[:, 2] = repeat(collect(1 : T), outer = length(m_names))
  
  output[(T * 0 + 1):(T * 1), 3:end] = get_pop_size(RR_ls, dg)
  output[(T * 1 + 1):(T * 2), 3:end] = get_pop_size(Rr_ls, dg)
  output[(T * 2 + 1):(T * 3), 3:end] = get_pop_size(rr_ls, dg)
  output[(T * 3 + 1):(T * 4), 3:end] = get_pro_R(RR_ls, Rr_ls, rr_ls)
  output[(T * 4 + 1):(T * 5), 3:end] = get_sur_rr(rr_ls, g_vals, s0, herb_ef, 
    g_pro, dg)
  output[(T * 5 + 1):(T * 6), 3:end] = get_TSR_adv(RR_ls, Rr_ls, rr_ls, dg,
    s0, herb_ef, g_pro, fec0, fec_cost, g_vals)
    
  return output
  
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

# produces a list that summaries the output_tup at each time speriod for %R and sur_rr
function trans_ex_summary(output_tup::Tuple{Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}},
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
  
  # get population in each cell at each time step, at each location for each G
  
  return snapshot
  
end

# function to produce a sinlge plot of %R over time and survival under mean g
function get_pro_R(RR::Array{Float64, 3}, Rr::Array{Float64, 3}, rr::Array{Float64, 3})
  
  num_RR = sum(RR, 1)
  num_RR = num_RR[1, :, :] # convert from 3D to 2D, note I actually want the transpose 
  num_Rr = sum(Rr, 1)
  num_Rr = num_Rr[1, :, :] 
  num_rr = sum(rr, 1)
  num_rr = num_rr[1, :, :] 
  
  return transpose((2 * num_RR + num_Rr) ./ (2 * (num_RR + num_Rr + num_rr)))
  
end

function get_pro_R(RR_ts::Array{Float64, 1}, Rr_ts::Array{Float64, 1}, 
  rr_ts::Array{Float64, 1})
  
  return (2 * RR_ts + Rr_ts) ./ (2 * (RR_ts + Rr_ts + rr_ts))
  
end

function get_pro_R(RR::Float64, Rr::Float64, rr::Float64)
  
  return (2 * RR + Rr) / (2 * (RR + Rr + rr))
  
end

# get population size
function get_pop_size(pop::Array{Float64, 3}, dg::Float64)

  tot_pop = sum(pop, 1) * dg
  
  return transpose(tot_pop[1, :, :]) # redice to 2D and transpose

end

# survival of TS susceptable population
function get_sur_rr(rr::Array{Float64, 3}, g_vals::Array{Float64, 1}, s0::Float64,
  herb_ef::Float64, g_pro::Float64, dg::Float64)

  pre_herb = get_pop_size(rr, dg)
  
  out = zeros(size(pre_herb))
  
  sur = 1 ./ (1 + exp(-(s0 - (herb_ef - min(herb_ef, g_vals * g_pro)))))

  sur_num = zeros(size(rr))
  
  T = size(rr)[3]
  X = size(rr)[2]
  
  for t in 1:T
  
    sur_num[:, :, t] = rr[:, :, t] .* sur
  
  end

  post_herb = get_pop_size(sur_num, dg)
  
  for t in 1:T
    for x in 1:X
  
      if pre_herb[t, x] > 0.0
	
	out[t, x] = post_herb[t, x] / pre_herb[t, x]
	
      end
      
    end
  end
  
  return out
  
end
# get fiotness of a given population, distribuited over g and space
function get_fitness(pop::Array{Float64, 3}, dg::Float64, sur_vec::Array{Float64, 1}, 
  fec_vec::Array{Float64, 1}, T::Int64)

  fit = zeros(size(pop))
 
  for t in 1:T
    
    fit[:, :, t] = pop[:, :, t] .* sur_vec .* fec_vec
 
  end
  
  fit_space_time = sum(fit, 1) * dg
  
  return transpose(fit_space_time[1, :, :])# reduce to 2D and transpose
  
end
# version for TSR 
function get_fitness(pop::Array{Float64, 3}, dg::Float64, sur::Float64, 
  fec_vec::Array{Float64, 1}, T::Int64)

  fit = zeros(size(pop))
 
  for t in 1:T
    
    fit[:, :, t] = pop[:, :, t] * sur .* fec_vec
 
  end
  
  fit_space_time = sum(fit, 1) * dg
  
  return transpose(fit_space_time[1, :, :])# redice to 2D and transpose
  
end

# get the realtive fitness advantage of TSR indivudals
function get_TSR_adv(RR_pop::Array{Float64, 3}, Rr_pop::Array{Float64, 3}, rr_pop::Array{Float64, 3}, dg::Float64,
  base_sur::Float64, h_eff::Float64, g_pro::Float64, fec0::Float64, fec_cost::Float64, g_vals::Array{Float64, 1})

  T = size(rr_pop)[3]
 
  g_effect_fec = 1 ./ (1 + exp(-(fec0 - abs(g_vals) * fec_cost)))
  sur_rr = 1 ./ (1 + exp(-(base_sur - (h_eff - min(h_eff, g_vals * g_pro)))))
  sur_R = 1 / (1 + exp(-(base_sur)))
  
  num_TSR = get_pop_size(RR_pop, dg) + get_pop_size(Rr_pop, dg)
  num_rr = get_pop_size(rr_pop, dg)
  
  KR = (get_fitness(RR_pop, dg, sur_R, g_effect_fec, T) + get_fitness(Rr_pop, dg, sur_R, g_effect_fec, T)) ./ num_TSR
  Krr = get_fitness(rr_pop, dg, sur_rr, g_effect_fec, T) ./ num_rr
  
  return KR ./ Krr
  
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
