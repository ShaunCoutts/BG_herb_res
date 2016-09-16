# simulation experiments for the TSR and NTSR joint evolution
# to test how TSR can invade into a population with NTSR

# Ideally I will be be able to specify which parameters to vary and
# the values to test them at. I could just pass every parameter 
# then I get length and run a for loop over it.

# int_?r_n should all be the same length. Each element is the number of black
# grass at each location. int_rr_1 passed unnamed as it will always be full.

function transplant_exper(int_rr_1::Array{Float64, 1}; int_rr_2::Array{Float64, 1}, 
  int_Rr_1::Array{Float64, 1} = [1, 2], int_Rr_2::Array{Float64, 1} = [1, 2], 
  herb_app_1:Array{Float64, 1} = [1, 2], herb_app_2::Array{Float64, 1} = [1,2],
  dx::Float64 = 1.0, dg::Float64 = 1.0, int_sd::Float64 = 1.4142, num_iter::Int64 = 50, 
  lower_g:Float64 = -10.0, upper_g::Float64 = 10.0, offspring_sd:Float64 = 1.0, 
  germ_prob = 0.5, fec_0 = 10.0, fec_cost = 10.0, fec_max = 60.0, dd_fec = 0.1, 
  base_sur = 10.0, herb_effect = 20.0, g_prot = 10.0, seed_sur = 0.5, pro_exposed = 0.8,
  scale_pollen = 32.0, shape_pollen = 3.32, seed_pro_short = 0.48, seed_mean_short = 0.58,
  seeds_to_mean_short = 0.44, seed_mean_long = 1.65, seeds_to_mean_long = 0.39)

  
  # TODO: the rest of the code
  
  
end