# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

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
function new_seeds_at_t!(RR_newseed::Array{Float64, 2}, Rr_newseed::Array{Float64, 2},
  rr_newseed::Array{Float64, 2},
  RR_mat::Array{Float64, 2}, Rr_mat::Array{Float64, 2}, rr_mat::Array{Float64, 2},
  RR_pollen::Array{Float64, 2}, Rr_pollen::Array{Float64, 2}, rr_pollen::Array{Float64, 2},
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1},
  fec_max = 100.0, dd_fec = 0.004, dg = 1.0)
  
  #holding array for density of new seeds new seeds before division
  RR_seeds = zeros(length(g_effect_fec))
  Rr_seeds = zeros(length(g_effect_fec))
  rr_seeds = zeros(length(g_effect_fec))
  new_seeds = zeros(length(g_effect_fec)) #generic holder vector to hold total seeds when they get split between G
  
  #divid fec_max by 3 since each maternal type will produce 3 seeds for each one that should exist, see why draw out all the combinations
  fec_corrected = fec_max / 3
  # hard code the G corsses, there is only 9 and and they will have fixed proportions: 
  for x in 1:size(RR_newseed)[2]
    
    num_at_x = (sum(RR_mat[:, x]) + sum(Rr_mat[:, x]) + sum(rr_mat[:, x])) * dg
    
    seeds_at_x = fec_corrected ./ (1 + g_effect_fec + dd_fec * num_at_x + dd_fec * num_at_x * g_effect_fec)
    #calcu numner of seeds for each G
    rr_seeds[:] = rr_mat[:, x] .* seeds_at_x
    Rr_seeds[:] = Rr_mat[:, x] .* seeds_at_x
    RR_seeds[:] = RR_mat[:, x] .* seeds_at_x
    
    # RR x RR seeds only produce RR seeds    
    RR_newseed[:, x] = g_mixing_kernel * 
      (RR_seeds[g_mixing_index[:, 1]] .* RR_pollen[g_mixing_index[:, 2], x]) * dg * dg 
    
    # RR x Rr seeds produce half RR seeds and half Rr seeds, number of seeds depends on maternal distrbution 
    new_seeds[:] = g_mixing_kernel * (RR_seeds[g_mixing_index[:, 1]] .* Rr_pollen[g_mixing_index[:, 2], x]) * dg * dg
    RR_newseed[:, x] = RR_newseed[:, x] + new_seeds * 0.5 #half the seeds produced are RR
    Rr_newseed[:, x] = new_seeds * 0.5 #other half go to Rr
    
    #Rr x RR
    new_seeds[:] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1]] .* RR_pollen[g_mixing_index[:, 2], x]) * dg * dg
    RR_newseed[:, x] = RR_newseed[:, x] + new_seeds * 0.5 #half the seeds produced are RR
    Rr_newseed[:, x] = Rr_newseed[:, x] + new_seeds * 0.5 #other half go to Rr
    
    #RR x rr produces only seeds of Rr
    Rr_newseed[:, x] = Rr_newseed[:, x] + g_mixing_kernel * 
      (RR_seeds[g_mixing_index[:, 1]] .* rr_pollen[g_mixing_index[:, 2], x]) * dg * dg
    #rr x RR
    Rr_newseed[:, x] = Rr_newseed[:, x] + g_mixing_kernel * 
      (rr_seeds[g_mixing_index[:, 1]] .* RR_pollen[g_mixing_index[:, 2], x]) * dg * dg
      
    #Rr x Rr produces all three genotypes  
    new_seeds[:] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1]] .* Rr_pollen[g_mixing_index[:, 2], x]) * dg * dg
    RR_newseed[:, x] = RR_newseed[:, x] + new_seeds * 0.25
    Rr_newseed[:, x] = Rr_newseed[:, x] + new_seeds * 0.5
    rr_newseed[:, x] = new_seeds * 0.25
    
    #Rr x rr produces Rr and rr seeds
    new_seeds[:] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1]] .* rr_pollen[g_mixing_index[:, 2], x]) * dg * dg
    rr_newseed[:, x] = rr_newseed[:, x] + new_seeds * 0.5
    Rr_newseed[:, x] = Rr_newseed[:, x] + new_seeds * 0.5
    #rr x Rr
    new_seeds[:] = g_mixing_kernel * (rr_seeds[g_mixing_index[:, 1]] .* Rr_pollen[g_mixing_index[:, 2], x]) * dg * dg
    rr_newseed[:, x] = rr_newseed[:, x] + new_seeds * 0.5
    Rr_newseed[:, x] = Rr_newseed[:, x] + new_seeds * 0.5
   
    #rr x rr produces only rr
    rr_newseed[:, x] = rr_newseed[:, x] + g_mixing_kernel * 
      (rr_seeds[g_mixing_index[:, 1]] .* rr_pollen[g_mixing_index[:, 2], x]) * dg * dg
  end

  return nothing
  
end

#TODO: test a matrix multiplication version of this 
function new_seeds_at_t_mm!(RR_newseed::Array{Float64, 2}, Rr_newseed::Array{Float64, 2},
  rr_newseed::Array{Float64, 2},
  RR_mat::Array{Float64, 2}, Rr_mat::Array{Float64, 2}, rr_mat::Array{Float64, 2},
  RR_pollen::Array{Float64, 2}, Rr_pollen::Array{Float64, 2}, rr_pollen::Array{Float64, 2},
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1},
  fec_max = 100.0, dd_fec = 0.004, dg = 1.0)
  
  #holding array for density of new seeds new seeds before division
  RR_seeds = zeros(size(RR_newseed)[1], size(RR_newseed)[2])
  Rr_seeds = zeros(size(Rr_newseed)[1], size(Rr_newseed)[2])
  rr_seeds = zeros(size(rr_newseed)[1], size(rr_newseed)[2])
  new_seeds = zeros(size(RR_newseed)[1], size(RR_newseed)[2]) #generic holder vector to hold total seeds when they get split between G
  
  #divid fec_max by 3 since each maternal type will produce 3 seeds for each one that should exist, see why draw out all the combinations
  fec_corrected = fec_max / 3
  # hard code the G corsses, there is only 9 and and they will have fixed proportions: 
    
  num_at_x = (sum(RR_mat, 1) + sum(Rr_mat, 1) + sum(rr_mat, 1)) * dg
    
  seeds_at_x = fec_corrected ./ (1 + g_effect_fec .+ dd_fec * num_at_x .+ dd_fec * (g_effect_fec * num_at_x))	
  
  #calcu numner of seeds for each G
  rr_seeds[:, :] = rr_mat .* seeds_at_x
  Rr_seeds[:, :] = Rr_mat .* seeds_at_x
  RR_seeds[:, :] = RR_mat .* seeds_at_x
    
  # RR x RR seeds only produce RR seeds    
  RR_newseed[:, :] = g_mixing_kernel * 
    (RR_seeds[g_mixing_index[:, 1], :] .* RR_pollen[g_mixing_index[:, 2], :]) * dg * dg 
    
  # RR x Rr seeds produce half RR seeds and half Rr seeds, number of seeds depends on maternal distrbution 
  new_seeds[:, :] = g_mixing_kernel * (RR_seeds[g_mixing_index[:, 1], :] .* Rr_pollen[g_mixing_index[:, 2], :]) * dg * dg
  RR_newseed[:, :] = RR_newseed + new_seeds * 0.5 #half the seeds produced are RR
  Rr_newseed[:, :] = new_seeds * 0.5 #other half go to Rr
    
  #Rr x RR
  new_seeds[:, :] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1], :] .* RR_pollen[g_mixing_index[:, 2], :]) * dg * dg
  RR_newseed[:, :] = RR_newseed + new_seeds * 0.5 #half the seeds produced are RR
  Rr_newseed[:, :] = Rr_newseed + new_seeds * 0.5 #other half go to Rr
    
  #RR x rr produces only seeds of Rr
  Rr_newseed[:, :] = Rr_newseed + g_mixing_kernel * 
    (RR_seeds[g_mixing_index[:, 1], :] .* rr_pollen[g_mixing_index[:, 2], :]) * dg * dg
  #rr x RR
  Rr_newseed[:, :] = Rr_newseed + g_mixing_kernel * 
    (rr_seeds[g_mixing_index[:, 1], :] .* RR_pollen[g_mixing_index[:, 2], :]) * dg * dg
    
  #Rr x Rr produces all three genotypes  
  new_seeds[:, :] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1], :] .* Rr_pollen[g_mixing_index[:, 2], :]) * dg * dg
  RR_newseed[:, :] = RR_newseed + new_seeds * 0.25
  Rr_newseed[:, :] = Rr_newseed + new_seeds * 0.5
  rr_newseed[:, :] = new_seeds * 0.25
  
  #Rr x rr produces Rr and rr seeds
  new_seeds[:, :] = g_mixing_kernel * (Rr_seeds[g_mixing_index[:, 1], :] .* rr_pollen[g_mixing_index[:, 2], :]) * dg * dg
  rr_newseed[:, :] = rr_newseed + new_seeds * 0.5
  Rr_newseed[:, :] = Rr_newseed + new_seeds * 0.5
  #rr x Rr
  new_seeds[:, :] = g_mixing_kernel * (rr_seeds[g_mixing_index[:, 1], :] .* Rr_pollen[g_mixing_index[:, 2], :]) * dg * dg
  rr_newseed[:, :] = rr_newseed + new_seeds * 0.5
  Rr_newseed[:, :] = Rr_newseed + new_seeds * 0.5
  
  #rr x rr produces only rr
  rr_newseed[:, :] = rr_newseed + g_mixing_kernel * 
    (rr_seeds[g_mixing_index[:, 1], :] .* rr_pollen[g_mixing_index[:, 2], :]) * dg * dg

  return nothing
  
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
function survival_at_t!(pop_at_t::Array{Float64, 2}, resist_G::Array{ASCIIString, 1}, G::ASCIIString, 
  herb_application::Array{Int64, 1}, sur_tup::Tuple{Float64, Array{Float64, 1}})
  
  if G in resist_G
      
    pop_at_t[:, :] = pop_at_t * sur_tup[1] 
    
  else
    
    for x in 1:size(pop_at_t)[2]
      pop_at_t[:, x] = pop_at_t[:, x] .* sur_tup[herb_application[x]] 
    end 
    
  end
  
  return nothing
  
end

# moves seeds to the next time step, kills the seeds in the process 
function seedbank_update!(seedbank_next::Array{Float64, 2}, seedbank_now::Array{Float64, 2}, seed_sur::Float64)

  seedbank_next[:, :] = seedbank_now * seed_sur
  
end

# Creates new plant and removes the germinated seeds from the seed bank
function new_plants!(ag_plants::Array{Float64, 2}, seed_bank::Array{Float64, 2}, germ_prob::Float64)
  
  ag_plants[:, :] = seed_bank * germ_prob
  seed_bank[:, :] = seed_bank * (1 - germ_prob)
  
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


