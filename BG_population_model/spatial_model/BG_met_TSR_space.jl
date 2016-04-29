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
  end
  
  return nothing
  
end

# Survival function for the whole landscape at a given time step t
# uses in place mutation. Note herb_application should be the same length as size(pop_at_t)[2] 
# and and g_vals should be the same length as size(pop_at_t)[1]
function survival_at_t!(pop_at_t::Array{Float64, 2}, base_sur::Float64, g_vals::Array{Float64, 1},
  resist_G::Array{ASCIIString, 1}, G::ASCIIString, herb_application::Array{Float64, 1}, herb_effect::Float64, g_prot::Float64, 
  pro_exposed::Float64)
  
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
# (intial_pop_size::Float32, intial_pop_g::Float32, int_pop_x::Array{Int8, 1})
function multi_iter(int_pop_RR::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]),
  int_pop_Rr::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]),
  int_pop_rr::Tuple{Float64, Float64, Array{Int64}, 1} = (0, 0, [4, 5, 6]); 
  int_sd::Float64 = 1.41, num_iter = 10, landscape_size = 10, space_res = 1, g_res = 1, lower_g = -10, 
  upper_g = 10, seed_expand = 3, germ_prob = 0.7)

  seed_disp_mat_2D = zeros(convert(Int32, ((landscape_size / space_res) + 1) ^ 2), 
    convert(Int32, ((landscape_size / space_res) + 1) ^ 2))
    
  seed_disp_mat_1D = zeros((landscape_size / space_res) + 1, (landscape_size / space_res) + 1)
    
  #build the above ground and seedbank evaluation points 
  g_ag_sb = g_points_builder(lower_g, upper_g, g_res, seed_expand)
 
  #set aside a chunck of memory for the landscapes for each genotype 
  RR_landscape = Array{Array{Float64, 2}}(num_iter)
  Rr_landscape = Array{Array{Float64, 2}}(num_iter)
  rr_landscape = Array{Array{Float64, 2}}(num_iter)
  # build a landscape for each genotype, for each time step 
  # these are 2D arrays of landscape_size x (upper_g - lower_g) * seed_expand 
  for t in 1:num_iter
    RR_landscape[t] = zeros(length(g_ag_sb[2]), landscape_size)
    Rr_landscape[t] = zeros(length(g_ag_sb[2]), landscape_size)
    rr_landscape[t] = zeros(length(g_ag_sb[2]), landscape_size)
  end
    
  #set the intial populations
  int_RR_dist = Normal(int_pop_RR[2], int_sd)
  int_Rr_dist = Normal(int_pop_Rr[2], int_sd)
  int_rr_dist = Normal(int_pop_rr[2], int_sd)
  #sets the intial population for each G at the locations specified in int_pop_x_G. 
  for x in int_pop_RR[3]
    RR_landscape[1][:, x] = pdf(int_RR_dist, g_ag_sb[2]) * int_pop_RR[1]
  end
  for x in int_pop_Rr[3]
    Rr_landscape[1][:, x] = pdf(int_Rr_dist, g_ag_sb[2]) * int_pop_Rr[1]
  end
  for x in int_pop_rr[3]
    rr_landscape[1][:, x] = pdf(int_rr_dist, g_ag_sb[2]) * int_pop_rr[1]
  end
  
  #a set of matrices to hold the above ground populations
  RR_ab_pop = zeros(length(g_ag_sb[1]), landscape_size)
  Rr_ab_pop = zeros(length(g_ag_sb[1]), landscape_size)
  rr_ab_pop = zeros(length(g_ag_sb[1]), landscape_size)
 
  # iterate through the timesteps
  for t in 2:num_iter
    
    new_plants!(RR_ab_pop, RR_landscape[g_ag_sb[3], :], germ_prob)
  
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
  
  
  

end 