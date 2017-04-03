# managment cost and effect models, also the immediate reward functions for econ only and econ + social
# will prodive arguments, and recive output from 
function state_update(state::Tuple{{Float64, 1}, Float64, Float64, Float64}, pop_sd::Float64,
  g_vals::Array{Float64, 1}, g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, 
  g_effect_fec::Array{Float64, 1}, seed_sur::Float64, germ::Float64, fec_max_corr::Float64, 
  sur_tup::Tuple{Float64, Array{Float64, 1}}, resist_G::Array{String, 1}, G::String, herb_application::Int64, spot_con::Bool, spot_mort::Float64, 
  dd_fec::Float64, dg::Float64)

# take an action and a state, return the next state and the immediate reward of going from s to s'
# action[1] = herb, action[2] = spot, action[3] = crop, action[4] = seedbank
function state_action2reward(state::Tuple{{Float64, 1}, Float64, Float64, Float64}, 
  action::Array{Int64, 1},  
  action_effect::Tuple{Array{Int64, 1}, Array{Bool, 1}, Array{Int64, 1}, Array{Float64, 1}},
  action_cost::Tuple{Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}},
  sur_tup::Tuple{Tuple{Float64, Array{Float64, 1}}, Tuple{Float64, Array{Float64, 1}}},  
  pop_sd::Float64, g_vals::Array{Float64, 1}, g_mixing_kernel::Array{Float64, 2}, spot_mort::Float64,
  g_mixing_index::Array{Int64, 2}, g_effect_fec::Array{Float64, 1}, fec_max_corr::Float64,
  resist_G::Array{String, 1}, dd_fec::Float64, dg::Float64)
  
  
  
  # state update
  if(action[CROP] != FAL) # for winter wheat or alternative crop
  
    new_pop = state_update(state, pop_sd, g_vals, g_mixing_kernel, g_mixing_index, 
      g_effect_fec, action_effect[BANK][action[BANK]], germ, fec_max_corr, 
      sur_tup[CROP], G, action_effect[HERB][action[HERB]], action_effect[SPOT][action[SPOT]], 
      spot_mort, dd_fec, dg)
  
  
  else # for fallow
  
    new_pop = state_update_fallow(state, pop_sd, g_vals, action_effect[BANK][action[BANK]], germ)
  
  end
  
  
  # these come from action state

end



# Builds the action space 
function build_action_space()

  # first build the wheat case 
  wheat_acts = hcat(repeat([NO_H, ONE_H, MUL_H], inner = 4),
    repeat(repeat([NO_SP, YS_SP], inner = 2), outer = 3),
    repeat([WHEAT], outer = 12),
    repeat([NO_SB, YS_SB], outer = 6))
    
  # built the ALT case
  alt_acts = hcat(repeat([NO_H, ONE_H, MUL_H], inner = 4),
    repeat(repeat([NO_SP, YS_SP], inner = 2), outer = 3),
    repeat([ALT], outer = 12),
    repeat([NO_SB, YS_SB], outer = 6))
    
  # Fallow is a special case becuase no point in herbicide or spot spray
  fal_acts = hcat(repeat([NO_H], outer = 2),
    repeat([NO_SP], outer = 2),
    repeat([FAL], outer = 2),
    [NO_SB, YS_SB])
  
  return vcat(wheat_acts, alt_acts, fal_acts)
  
end