# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

using Distributions

# tatkes a lower and upper domain limit and resolution and expantion factor for the seedbank 
# returns a tuple of three other tuples
# 1: above ground values of g being evaluated 
# 2: seedabank values of g being evaluated
# 3: indexes of above ground values of g in the tuple of seedbank values of g
function g_points_builder(low_point, upper_point, res, seed_expand)
  above_ground_eval = (low_point : res : upper_point)
  seed_eval = (low_point * seed_expand : res : upper_point * seed_expand)
  above_ground_ind = tuple(findin(seed_eval, above_ground_eval)...)
  
  return (above_ground_eval, seed_eval, above_ground_ind)
end

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
function new_seeds_at_x!(RR_newseed::Array{Float64, 2}, Rr_newseed::Array{Float64, 2}
  rr_newseed::Array{Float64, 2},
  RR_mat::Array{Float64, 2}, Rr_mat::Array{Float64, 2}, rr_mat::Array{Float64, 2},
  RR_pollen::Array{Float64, 2}, Rr_pollen::Array{Float64, 2}, rr_pollen::Array{Float64, 2},
  g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1};
  fec_max = 100.0, dd_fec = 0.004, dg = 0.5)
  
  #holding array for density of new seeds new seeds before division
  RR_seeds = zeros(length(g_effect_fec))
  Rr_seeds = zeros(length(g_effect_fec))
  rr_seeds = zeros(length(g_effect_fec))
  new_seeds = zeros(length(g_effect_fec)) #generic holder vector to hold total seeds when they get split between G
  
  #divid fec_max by 3 since each maternal type will produce 3 seeds for each one that should exist, see why draw out all the combinations
  fec_corrected = fec_max / 3
  # hard code the G corsses, there is only 9 and and they will have fixed proportions: 
  for x in 1:size(RR_newseed)[2]
    
    num_at_x = (sum(RR_mat) + sum(Rr_mat) + sum(rr_mat)) * dg
    
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

# disperse the seeds produced at each location to every other location, suming the number of seeds of 
# each G that arrive at each location 






# Survival function that takes an element of population array and reduces it in place
# be aware that pop_at_x and g_vals need to match up, that s one element in pop_at_x
# should corospond to a element of g_vals 
function survival_at_x!(pop_at_x::Array{Float64, 1}, base_sur::Float64, g_vals::Array{Float64, 1} 
  resist_G::Array{ASCIIString, 1}, G::ASCIIString, herb_application::Float64, herb_effect::Float64, g_prot::Float64, 
  pro_exposed::Float64)
  
  if G in resist_G
    pop_at_x = pop_at_x * (1 / (1 + exp(-base_sur)))
  else
    pop_at_x = pop_at_x .* ((1 - pro_exposed) / (1 + exp(-base_sur))) + 
      (pro_exposed ./ (1 + exp(-base_sur - herb_application * 
      (herb_effect - min(herb_effect, g_vals * g_prot)))))
  end`
  
  return nothing
  
end

# Survival function for the whole landscape at a given time step t
# uses in place mutation. Note herb_application should be the same length as size(pop_at_t)[2] 
# and and g_vals should be the same length as size(pop_at_t)[1]
function survival_at_t!(pop_at_t::Array{Float64, 2}, base_sur::Float64, g_vals::Array{Float64, 1},
  resist_G::Array{ASCIIString, 1}, G::ASCIIString, herb_application::Array{Float64, 1}, herb_effect::Float64, 
  g_prot::Float64, pro_exposed::Float64)
  
  for x in 1:size(pop_at_t)[2]
    survival_at_x(pop_at_t[:, x], base_sur, g_vals, resist_G, G, herb_application[x],
      herb_effect, g_prot, pro_exposed)
  end
  
  return nothing
  
end

function new_plants!(ag_plants::Array{Float64, 2}, seed_bank::Array{Float64, 2}, germ_prob)
  
  ag_plants[:, :] = seed_bank[:, :] * germ_prob
  
  return nothing
end

# Takes the current state of the population -> next state of the population  
function single_iter(current_pop::Array{Float64, 2}, next_pop::Array{Float64, 2}; seed_disp_mat::Array{Float64, 2}, 
  pollen_disp_mat::Array{Float64, 2})
 
 
 return next_pop 
end

# the first three positional arguments (int_pop_RR, int_pop_Rr, int_pop_rr) are tuples of
# (intial_pop_size::Float64, intial_pop_g::Float64, int_pop_x::Array{Int8, 1})
function multi_iter_1D(int_pop_RR::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]),
  int_pop_Rr::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]),
  int_pop_rr::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]); 
  int_sd::Float64 = 1.41, num_iter = 10, landscape_size = 10, space_res = 1, g_res = 1, lower_g = -10, 
  upper_g = 10, seed_expand = 2, germ_prob = 0.7, offspring_sd = 1.0, fec0 = 10, fec_cost = 2)

    
  seed_disp_mat_1D = zeros((landscape_size / space_res) + 1, (landscape_size / space_res) + 1)
    
  # build the metabolic resistance evaluation points 
  g_vals = (lower_g : g_res : upper_g)
 
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
  g_effect_fec = exp(-(fec0 - g_vals * fec_cost))
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
  
  # iterate through the timesteps
  for t in 2:num_iter
    
    new_plants!(RR_ab_pop, RR_landscape[t][g_vals, :], germ_prob)
    new_plants!(Rr_ab_pop, Rr_landscape[t][g_vals, :], germ_prob)
    new_plants!(rr_ab_pop, rr_landscape[t][g_vals, :], germ_prob)
    
    ## pollen at arrived each location for each g from all other locations  
    pollen_RR[:, :] = RR_ab_pop * pollen_disp_mat
    pollen_Rr[:, :] = Rr_ab_pop * pollen_disp_mat
    pollen_rr[:, :] = rr_ab_pop * pollen_disp_mat
    total_pollen[:] = sum(pollen_RR, 1) + sum(pollen_Rr, 1) + sum(pollen_rr, 1)
    #normalise the pollen counts, double intergration across g and and all x as pollen comes from all x
    pollen_RR[:, :] = pollen_RR ./ (total_pollen * g_res * space_res)
    pollen_Rr[:, :] = pollen_Rr ./ (total_pollen * g_res * space_res)
    pollen_rr[:, :] = pollen_rr ./ (total_pollen * g_res * space_res)
    
    
    #create new seeds
    
    
    # disperse the seeds
    RR_landscape[t][:, :] = RR_landscape  + (RR_newseed * seed_disp_mat_1D) * space_res 
    
    
    ## TODO test I need the double intergration in the pollen counts I think I do
    
    
    
    
    
  end
  
 
  
  @time seed_disp_mat_builder_1D!(seed_disp_mat, res = space_res)
  
  pollen_disp_mat = zeros(convert(Int32, (landscape_size / space_res) + 1), 
    convert(Int32, (landscape_size / space_res) + 1))
  
  @time pollen_disp_mat_builder_1D!(pollen_disp_mat, res = space_res)
  
  
  @time seed_disp_mat_builder_2D!(seed_disp_mat_2D, res = space_res, pro_short = 0.48, mean_dist_short = 0.58, 
    pro_seeds_to_mean_short = 0.44, mean_dist_long = 1.65, pro_seeds_to_mean_long = 0.39)
  


end


# function to run and call the the other functions and scripts, will eventually run the whole thing
function main_calling_function()
  cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model")
  include("BG_met_TSR_space_dispersal_functions.jl")
  srand(321) #set random seed
  
  
#   seed_disp_mat_2D = zeros(convert(Int32, ((landscape_size / space_res) + 1) ^ 2), 
#     convert(Int32, ((landscape_size / space_res) + 1) ^ 2))
#   

end 