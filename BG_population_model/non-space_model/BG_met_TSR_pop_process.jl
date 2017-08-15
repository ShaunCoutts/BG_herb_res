# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

using Distributions

# fills the g_mixing kernel with a normal offspring dist for each coloumn
# witha mean of g_m * g_p for every combination
function fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)

  m_p_comb = 1
  for g_m in g_vals
    for g_p in g_vals
      g_mixing_kernel[:, m_p_comb] = pdf(Normal(0.5 * g_m + 0.5 * g_p, offspring_sd), g_vals)
      
      m_p_comb += 1
    end
  end

end

#taks a vector and resturs every pair of indicies
function generate_index_pairs(vect)
  inds = Array(Int64, length(vect) ^ 2, 2)
  count = 1
  for i in eachindex(vect)
    for j in eachindex(vect)
      inds[count, :] = [i j]
      count += 1
    end
  end
  
  return inds
  
end

# mixing of g and G at a given location
# fist three argumaents are the arrays to save the new RR, Rr and rr seed in for a given location
# mext three are the maternal distirbution i.e. survivours at x, for each G
# next three is the pollen dist that arrived at x
# followed by the mixing kernel for g, and the indicies to apply to the parent dists for that mixing 
function new_seeds_at_t!(RR_newseed::Array{Float64, 1}, Rr_newseed::Array{Float64, 1},
  rr_newseed::Array{Float64, 1},
  RR_mat::Array{Float64, 1}, Rr_mat::Array{Float64, 1}, rr_mat::Array{Float64, 1},
  RR_pollen::Array{Float64, 1}, Rr_pollen::Array{Float64, 1}, rr_pollen::Array{Float64, 1},
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, fec_max::Float64, dg::Float64)
  
  
  # hard code the G corsses, there is only 9 and and they will have fixed proportions: 
    
  #calcu numner of seeds for each G
  rr_seeds = rr_mat * fec_max
  Rr_seeds = Rr_mat * fec_max
  RR_seeds = RR_mat * fec_max
  new_seeds = zeros(size(RR_mat)[1]) #generic holder vector to hold total seeds when they get split between G
    
  # RR x RR seeds only produce RR seeds    
  RR_newseed[:] = g_mixing_kernel * 
    (RR_seeds[g_mixing_index[:, 1]] .* RR_pollen[g_mixing_index[:, 2]]) * dg * dg 
    
  # RR x Rr seeds produce half RR seeds and half Rr seeds, number of seeds depends on maternal distrbution 
  new_seeds[:] = g_mixing_kernel * (RR_seeds[g_mixing_index[:, 1]] .* Rr_pollen[g_mixing_index[:, 2]]) * dg * dg
  RR_newseed[:] = RR_newseed[:] + new_seeds * 0.5 #half the seeds produced are RR
  Rr_newseed[:] = new_seeds * 0.5 #other half go to Rr
    
  #Rr x RR
  new_seeds[:] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1]] .* RR_pollen[g_mixing_index[:, 2]]) * dg * dg
  RR_newseed[:] = RR_newseed[:] + new_seeds * 0.5 #half the seeds produced are RR
  Rr_newseed[:] = Rr_newseed[:] + new_seeds * 0.5 #other half go to Rr
  
  #RR x rr produces only seeds of Rr
  Rr_newseed[:] = Rr_newseed[:] + g_mixing_kernel * 
    (RR_seeds[g_mixing_index[:, 1]] .* rr_pollen[g_mixing_index[:, 2]]) * dg * dg
  #rr x RR
  Rr_newseed[:] = Rr_newseed[:] + g_mixing_kernel * 
    (rr_seeds[g_mixing_index[:, 1]] .* RR_pollen[g_mixing_index[:, 2]]) * dg * dg
    
  #Rr x Rr produces all three genotypes  
  new_seeds[:] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1]] .* Rr_pollen[g_mixing_index[:, 2]]) * dg * dg
  RR_newseed[:] = RR_newseed[:] + new_seeds * 0.25
  Rr_newseed[:] = Rr_newseed[:] + new_seeds * 0.5
  rr_newseed[:] = new_seeds * 0.25
  
  #Rr x rr produces Rr and rr seeds
  new_seeds[:] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1]] .* rr_pollen[g_mixing_index[:, 2]]) * dg * dg
  rr_newseed[:] = rr_newseed[:] + new_seeds * 0.5
  Rr_newseed[:] = Rr_newseed[:] + new_seeds * 0.5
  #rr x Rr
  new_seeds[:] = g_mixing_kernel * (rr_seeds[g_mixing_index[:, 1]] .* Rr_pollen[g_mixing_index[:, 2]]) * dg * dg
  rr_newseed[:] = rr_newseed[:] + new_seeds * 0.5
  Rr_newseed[:] = Rr_newseed[:] + new_seeds * 0.5
  
  #rr x rr produces only rr
  rr_newseed[:] = rr_newseed[:] + g_mixing_kernel * 
    (rr_seeds[g_mixing_index[:, 1]] .* rr_pollen[g_mixing_index[:, 2]]) * dg * dg

  return nothing
  
end

# density effect 
function density_effect(f_d::Float64, tot_pop::Float64)
  
  return 1 / (1 + f_d * tot_pop)

end

# pre-calc resistance cost 
function resist_cost_pre_calc(f0::Float64, f_cost::Float64, g_vals::Array{Float64, 1})

  return 1 ./ (1 + exp(-(f0 - f_cost * abs(g_vals)))) 

end

# calculate the two survival rates used, one a scalar for none exposed plants and 
# one for exposed plants
# be aware that pop_at_x and g_vals need to match up, that s one element in pop_at_x
# should corospond to a element of g_vals 
function survival_pre_calc(base_sur::Float64, g_vals::Array{Float64, 1}, herb_effect::Float64, 
  g_prot::Float64, pro_exposed::Float64)
  
  sur_non_exposed = 1 / (1 + exp(-base_sur))
  
  sur_exposed = ((1 - pro_exposed) / (1 + exp(-base_sur))) + 
      (pro_exposed ./ (1 + exp(-(base_sur - (herb_effect - min(herb_effect, g_vals * g_prot))))))
  
  return (sur_non_exposed, sur_exposed)
   
end

# Survival function for the whole landscape at a given time step t
# uses in place mutation. Note herb_application should be the same length as size(pop_at_t)[2] 
# and and g_vals should be the same length as size(pop_at_t)[1]
function survival_at_t!(pop_at_t::Array{Float64, 1}, resist_G::Array{String, 1}, G::String, 
  herb_application::Array{Int64, 1}, sur_tup::Tuple{Float64, Array{Float64, 1}})
  
  if G in resist_G
      
    pop_at_t[:] = pop_at_t[:] * sur_tup[1] 
    
  else
    
    pop_at_t[:] = pop_at_t[:] .* sur_tup[herb_application[x]] 
    
  end
  
  return nothing
  
end

# moves seeds to the next time step, kills the seeds in the process 
function seedbank_update!(seedbank_next::Array{Float64, 1}, seedbank_now::Array{Float64, 1}, seed_sur::Float64)

  seedbank_next[:] = seedbank_now * seed_sur
  
end

# Creates new plant and removes the germinated seeds from the seed bank
function new_plants!(ag_plants::Array{Float64, 1}, seed_bank::Array{Float64, 1}, germ_prob::Float64)
  
  ag_plants[:] = seed_bank * germ_prob
  seed_bank[:] = seed_bank * (1 - germ_prob)
  
  return nothing
  
end

#get the mean and sd of a location
function dist_summary(dist::Array{Float64, 1}, g_vals::Array{Float64, 1}, dg::Float64)
  
  total_sum = sum(dist) * dg
  if total_sum > 0.000001
    approx_mean = sum((dist / total_sum) .* g_vals) * dg
    approx_sd = sqrt(sum(((g_vals - approx_mean) .^ 2) .* (dist / total_sum)) * dg)
  else
    approx_mean = 0.0
    approx_sd = 1.0
  end
  
  return [approx_mean, approx_sd, total_sum]

end

# function to get the output of interest for the parameter filtering (total above 
# ground pop after herbicide) from the multi-iteration output (mean, SD and total 
# seed bank at each location at each time).
# INPUT: pop_snapshot is a 2d slice of the 3D array produced by multi_iter_1D()
# OUTPUT: vector of population metrics
function pop_summaries(pop_snapshot::Array{Float64, 1}, g_ext::Array{Float64, 1}, dg::Float64, param_var::Array{Float64, 1}, 
  base_sur::Float64)

  # extract some parameters that are useful
  germ_prob = param_var[4]
  herb_effect = param_var[9]
  g_prot = param_var[10]
  pro_exposed = param_var[12]
  
  # in order to calculate several of the population summaries I need to reconstruct the population
  g_vals = collect(g_ext[1] : dg : g_ext[2])
  pop_sb_RR = zeros(length(g_vals)) # blank landscape seedbank
  pop_sb_Rr = zeros(length(g_vals)) 
  pop_sb_rr = zeros(length(g_vals)) 
  
  #creat the seed bank and above ground numbers both pre and post herbicide  
  pop_sb_RR[:] = pdf(Normal(pop_snapshot[1], pop_snapshot[2]), g_vals) * pop_snapshot[3]
  pop_sb_Rr[:] = pdf(Normal(pop_snapshot[4], pop_snapshot[5]), g_vals) * pop_snapshot[6]
  pop_sb_rr[:] = pdf(Normal(pop_snapshot[7], pop_snapshot[8]), g_vals) * pop_snapshot[9]
 
  # create the above ground populations, the post herbcide matricies also creatged to be mutated later
  pop_ab_RR = pop_sb_RR * germ_prob 
  pop_ab_Rr = pop_sb_Rr * germ_prob
  pop_ab_rr = pop_sb_rr * germ_prob 
  
  pop_ab_ph_RR = deepcopy(pop_ab_RR)
  pop_ab_ph_Rr = deepcopy(pop_ab_Rr)
  pop_ab_ph_rr = deepcopy(pop_ab_rr)
  
  #do one iteration of survival to ge the number of post herbicide above ground plants
  sur_pre_calc = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)
  survival_at_t!(pop_ab_ph_RR, ["RR", "Rr"], "RR", convert(Array{Int64, 1}, ones(length(landscape)) + 1), sur_pre_calc)
  survival_at_t!(pop_ab_ph_Rr, ["RR", "Rr"], "Rr", convert(Array{Int64, 1}, ones(length(landscape)) + 1), sur_pre_calc)
  survival_at_t!(pop_ab_ph_rr, ["RR", "Rr"], "rr", convert(Array{Int64, 1}, ones(length(landscape)) + 1), sur_pre_calc)
  
  
  #population numbers
  num_sb_RR = sum(pop_snapshot[3, :]) * dx
  num_sb_Rr = sum(pop_snapshot[6, :]) * dx
  num_sb_rr = sum(pop_snapshot[9, :]) * dx
  num_sb_tot = num_sb_RR + num_sb_Rr + num_sb_rr
  
  num_ab_tot = (sum(pop_ab_RR) + sum(pop_ab_Rr) + sum(pop_ab_rr)) * dx * dg  
  num_ab_ph_tot = (sum(pop_ab_ph_RR) + sum(pop_ab_ph_Rr) + sum(pop_ab_ph_rr)) * dx * dg  
 
  #resistance numbers
  mean_g_RR = (sum(pop_snapshot[1, :] .* pop_snapshot[3, :]) * dx) / num_sb_RR 
  mean_g_Rr = (sum(pop_snapshot[4, :] .* pop_snapshot[6, :]) * dx) / num_sb_Rr 
  mean_g_rr = (sum(pop_snapshot[7, :] .* pop_snapshot[9, :]) * dx) / num_sb_rr 
  mean_g_pop = (mean_g_RR * num_sb_RR + mean_g_Rr * num_sb_Rr + mean_g_rr * num_sb_rr) / num_sb_tot
   
  # get the sperad indicies
  pro_RR_x = sum(sum(pop_sb_RR, 1) .> 1) / length(landscape)
  pro_Rr_x = sum(sum(pop_sb_Rr, 1) .> 1) / length(landscape)
  pro_rr_x = sum(sum(pop_sb_rr, 1) .> 1) / length(landscape) 
  pro_all_x = sum((sum(pop_sb_RR, 1) .> 1) | (sum(pop_sb_Rr, 1) .> 1) | (sum(pop_sb_rr, 1) .> 1)) / length(landscape)
  
  return [num_sb_RR, num_sb_Rr, num_sb_rr, num_sb_tot, num_ab_tot, num_ab_ph_tot, 
    mean_g_RR, mean_g_Rr, mean_g_rr, mean_g_pop, pro_RR_x, pro_Rr_x, pro_rr_x, pro_all_x]
 
end
