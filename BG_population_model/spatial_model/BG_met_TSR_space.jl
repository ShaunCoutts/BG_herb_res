# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

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
function survival_at_x!(pop_at_x::Array{Float64, 1}, base_sur, g_vals::Array{Float32, 1} 
  resist_G, G, herb_application, herb_effect, g_prot, pro_exposed)
  if G in resist_G
    pop_at_x = pop_at_x * (1 / (1 + exp(-base_sur)))
  else
    pop_at_x = pop_at_x .* ((1 - pro_exposed) / (1 + exp(-base_sur))) + 
      (pro_exposed ./ (1 + exp(-base_sur - herb_application * 
      (herb_effect - min(herb_effect, g_vals * g_prot)))))
  end
  
  return nothing
  
end



# Takes the current state of the population -> next state of the population  
function single_iter(current_pop::Array{Float64, 2}, next_pop::Array{Float64, 2}; seed_disp_mat::Array{Float64, 2}, 
  pollen_disp_mat::Array{Float64, 2})
 
 
 return next_pop 
end


function multi_iter(num_iter = 10, landscape_size = 10, space_res = 1, g_res = 1, lower_g = -10, 
  upper_g = 10, seed_expand = 3)

  seed_disp_mat_2D = zeros(convert(Int32, ((landscape_size / space_res) + 1) ^ 2), 
    convert(Int32, ((landscape_size / space_res) + 1) ^ 2))
    
  seed_disp_mat_1D = zeros((landscape_size / space_res) + 1, (landscape_size / space_res) + 1)
    
  #build the above ground and seedbank evaluation points 
  g_ag_sb = g_points_builder(lower_g, upper_g, res = g_res, seed_expand)
 
  #build a landscape for each genotype, these are 2D arrays of landscape_size x (upper_g - lower_g) * seed_expand 
  RR_landscape = zeros(landscape_size, length(a_ag_sb[2]))
  
  
  
  @time seed_disp_mat_builder_1D!(seed_disp_mat, res = space_res)
  
  pollen_disp_mat = zeros(convert(Int32, (landscape_size / space_res) + 1), 
    convert(Int32, (landscape_size / space_res) + 1))
  
  @time pollen_disp_mat_builder_1D!(pollen_disp_mat, res = space_res)
  
  
  @time seed_disp_mat_builder_2D!(seed_disp_mat, res = space_res, pro_short = 0.48, mean_dist_short = 0.58, 
    pro_seeds_to_mean_short = 0.44, mean_dist_long = 1.65, pro_seeds_to_mean_long = 0.39)
  


end


# function to run and call the the other functions and scripts, will eventually run the whole thing
function main_calling_function()
  cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model")
  include("BG_met_TSR_space_dispersal_functions.jl")
  
  
  

end 