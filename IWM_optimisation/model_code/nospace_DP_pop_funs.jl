# populaiton  model functions for dynamic programming 
# these should take a state (Array{Any, 1}), build a population from that 
# state, step the population foward one step, calcualte the value at that new state. 

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

# take a array of seed bank (distributed over g), update to new seedbank 
function sb_update!(RR_sb::Array{Float64, 1}, Rr_sb::Array{Float64, 1},
  rr_sb::Array{Float64, 1}, tot_newseed::Array{Float64, 1},
  RR_mat::Array{Float64, 1}, Rr_mat::Array{Float64, 1}, rr_mat::Array{Float64, 1},
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1},
  fec_max_corr::Float64, dd_fec::Float64, dg::Float64)
  
  # get num seeds per mother and also calculate pollen vectors
    
  tot_num = (sum(RR_mat) + sum(Rr_mat) + sum(rr_mat)) * dg
  
  seeds = fec_max_corr ./ (1.0 + g_effect_fec .+ (dd_fec * tot_num + dd_fec * (g_effect_fec * tot_num)))	
  
  RR_pollen = RR_mat[g_mixing_index[:, 2]] / tot_num
  Rr_pollen = Rr_mat[g_mixing_index[:, 2]] / tot_num
  rr_pollen = rr_mat[g_mixing_index[:, 2]] / tot_num
  
  #calcu numner of seeds for each G
  rr_seeds = rr_mat[g_mixing_index[:, 1]] .* seeds[g_mixing_index[:, 1]]
  Rr_seeds = Rr_mat[g_mixing_index[:, 1]] .* seeds[g_mixing_index[:, 1]]
  RR_seeds = RR_mat[g_mixing_index[:, 1]] .* seeds[g_mixing_index[:, 1]]
  #holding array for density of new seeds new seeds before division
  tot_newseed[:] = 0.0
  
  # g mixing works via the g_mixing kernel that is g_val by g_val^2 matrix, multiplied by the 
  # pollen and seed vectors, that are g_val^2 coloumn vectors, product is g_val by 1 matrix 
  # (i.e. coloumn vecrtor). relies on G_mat, G_sb and g_efect_fec all having g_val rows 
  
  # hard code the G corsses, there is only 9 and and they will have fixed proportions: 
  # RR x RR seeds only produce RR seeds    
  RR_sb[:] += g_mixing_kernel * (RR_seeds .* RR_pollen) * dg * dg 
    
  # RR x Rr seeds produce half RR seeds and half Rr seeds, number of seeds 
  # depends on maternal distrbution 
  tot_newseed[:] = g_mixing_kernel * (RR_seeds .* Rr_pollen) * dg * dg
  RR_sb[:] += tot_newseed * 0.5 #half the seeds produced are RR
  Rr_sb[:] += tot_newseed * 0.5 #other half go to Rr
    
  #Rr x RR
  tot_newseed[:] = g_mixing_kernel * (Rr_seeds .* RR_pollen) * dg * dg
  RR_sb[:] += tot_newseed * 0.5 #half the seeds produced are RR
  Rr_sb[:] += tot_newseed * 0.5 #other half go to Rr
    
  #RR x rr produces only seeds of Rr
  Rr_sb[:] += g_mixing_kernel * (RR_seeds .* rr_pollen) * dg * dg
  #rr x RR
  Rr_sb[:] += g_mixing_kernel * (rr_seeds .* RR_pollen) * dg * dg
    
  #Rr x Rr produces all three genotypes  
  tot_newseed[:] = g_mixing_kernel * (Rr_seeds .* Rr_pollen) * dg * dg
  RR_sb[:] += tot_newseed * 0.25
  Rr_sb[:] += tot_newseed * 0.5
  rr_sb[:] += tot_newseed * 0.25
  
  #Rr x rr produces Rr and rr seeds
  tot_newseed[:] = g_mixing_kernel * (Rr_seeds .* rr_pollen) * dg * dg
  rr_sb[:] += tot_newseed * 0.5
  Rr_sb[:] += tot_newseed * 0.5
  #rr x Rr
  tot_newseed[:] = g_mixing_kernel * (rr_seeds .* Rr_pollen) * dg * dg
  rr_sb[:] += tot_newseed * 0.5
  Rr_sb[:] += tot_newseed * 0.5
  
  #rr x rr produces only rr
  rr_sb[:] += g_mixing_kernel * (rr_seeds .* rr_pollen) * dg * dg

  return nothing
  
end

# Creates new plant and removes the germinated seeds from the seed bank
function new_plants(seed_bank::Array{Float64, 1}, germ_prob::Float64)
  
  ag_plants = seed_bank * germ_prob
  seed_bank[:] *= (1 - germ_prob)
  
  return ag_plants
end

# calculate the two survival rates used, one a scalar for none exposed plants and 
# one for exposed plants
# be aware that pop_at_x and g_vals need to match up, that s one element in pop
# should corospond to a element of g_vals 
function survival_pre_calc(base_surW::Float64, base_surA::Float64, 
  g_vals::Array{Float64, 1}, herb_effect::Float64, g_prot::Float64, pro_exposed::Float64)
  
  sur_non_exposedW = 1 / (1 + exp(-base_surW))
  sur_exposedW = (pro_exposed ./ (1 + exp(-(base_surW - (herb_effect - min(herb_effect, g_vals * g_prot))))))
  
  sur_non_exposedA = 1 / (1 + exp(-base_surA))
  sur_exposedA = (pro_exposed ./ (1 + exp(-(base_surA - (herb_effect - min(herb_effect, g_vals * g_prot))))))
  
  return ((sur_non_exposedW, sur_exposedW), (sur_non_exposedA, sur_exposedA))
   
end

# takes a population and applies herbicide survival to it
# herb_application determines how many times herbicide is 
# applied, can be 0 
function survival_herb!(pop::Array{Float64, 1}, G::String, 
  herb_application::Int64, sur_tup::Tuple{Float64, Array{Float64, 1}})
  
  #base sur the whole population is exposed to 
  pop[:] *= sur_tup[1]
  
  if !(G in resist_G)
    for i in 1:herb_application
      pop[:] = pop .* sur_tup[2] 
    end
  end
  
  return nothing
  
end

# Survival under spot spray where wheat dies too, nuclear option
function survival_spot!(pop::Array{Float64, 1}, spot::Bool, spot_mort::Float64)

  if(spot)
    pop[:] *= spot_mort
  end

  return nothing

end

# take a state (Tuple(4)), iterate one step foward to get the next population.
# state[1][1] = mean(RR_g); state[1][2] = mean(Rr_g); state[1][3] = mean(rr_g);
# state[2] = total_pop;, state[3] = pro_RR; state[4] = pro_Rr;
# make the assumption that pop_sd is constant across action, seedbank sizes
# return a 2D array, g_val by 3, one coloumn for each G, each holding the 
# seedbank (to get next state) and total above ground pop after all mortality
function state_update(state::Tuple{{Float64, 1}, Float64, Float64, Float64}, pop_sd::Float64,
  g_vals::Array{Float64, 1}, g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, 
  g_effect_fec::Array{Float64, 1}, seed_sur::Float64, germ::Float64, fec_max_corr::Float64, 
  sur_tup_act::Tuple{Float64, Array{Float64, 1}}, G::String, 
  herb_application::Int64, spot_con::Bool, spot_mort::Float64, dd_fec::Float64, dg::Float64)

  # unpack the state variable, 
  pro_RR = state[S_pRR]
  pro_Rr = (1 - state[S_pRR]) * state[S_pRr]
  pro_rr = 1 - pro_RR - pro_Rr
  # state[2] is total pop 
  RR_sb = pdf(Normal(state[S_g][S_gRR], pop_sd), g_vals) * pro_RR * state[S_N] 
  Rr_sb = pdf(Normal(state[S_g][S_gRr], pop_sd), g_vals) * pro_Rr * state[S_N] 
  rr_sb = pdf(Normal(state[S_g][S_grr], pop_sd), g_vals) * pro_rr * state[S_N] 
  # array to temp hold new seeds before they are dvided amoung G
  tot_newseed = zeros(size(RR_sb)[1])
  
  # seed survial
  RR_sb *= seed_sur
  Rr_sb *= seed_sur
  rr_sb *= seed_sur
    
  # germination
  RR_mat = new_plants(RR_sb, germ)
  Rr_mat = new_plants(Rr_sb, germ)
  rr_mat = new_plants(rr_sb, germ)
  
  # above ground survival 
  survival_herb!(RR_mat, resist_G, "RR", herb_application, sur_tup_act)
  survival_herb!(Rr_mat, resist_G, "Rr", herb_application, sur_tup_act)
  survival_herb!(rr_mat, resist_G, "rr", herb_application, sur_tup_act)
  
  tot_ab_post_herb = (sum(RR_mat) + sum(Rr_mat) + sum(rr_mat)) * dg
  
  # spot spray that kills wheat also 
  if(spot_con)
    RR_mat *= spot_mort
    Rr_mat *= spot_mort
    rr_mat *= spot_mort
  end
  
  tot_ab_post_spot = (sum(RR_mat) + sum(Rr_mat) + sum(rr_mat)) * dg
  
  sb_update!(RR_sb, Rr_sb, rr_sb, tot_newseed, RR_mat, Rr_mat, rr_mat,
    g_mixing_kernel, g_mixing_index, g_effect_fec, fec_max_corr, dd_fec, dg)
  
  return [RR_sb, Rr_sb, rr_sb, tot_ab_post_herb, tot_ab_post_spot] 
  
end

# State update when fallow action taken, we assume all above ground population
# killed after germination so only need seed survival and germination to get next state 
function state_update_fallow(state::Tuple{{Float64, 1}, Float64, Float64, Float64}, pop_sd::Float64,
  g_vals::Array{Float64, 1}, seed_sur::Float64, germ::Float64)

  # unpack the state variable, 
  pro_RR = state[S_pRR]
  pro_Rr = (1 - state[S_pRR]) * state[S_pRr]
  pro_rr = 1 - pro_RR - pro_Rr
  # state[2] is total pop 
  RR_sb = pdf(Normal(state[S_g][S_gRR], pop_sd), g_vals) * pro_RR * state[S_N] 
  Rr_sb = pdf(Normal(state[S_g][S_gRr], pop_sd), g_vals) * pro_Rr * state[S_N] 
  rr_sb = pdf(Normal(state[S_g][S_grr], pop_sd), g_vals) * pro_rr * state[S_N] 
  
  # seed survial
  RR_sb *= seed_sur
  Rr_sb *= seed_sur
  rr_sb *= seed_sur
    
  # germination
  RR_sb *= (1 - germ)
  Rr_sb *= (1 - germ)
  rr_sb *= (1 - germ)
  
  return [RR_sb, Rr_sb, rr_sb, 0.0, 0.0] 
  
end


# Takes output of state_update() and calculates the reward, use econ and econ +  social
