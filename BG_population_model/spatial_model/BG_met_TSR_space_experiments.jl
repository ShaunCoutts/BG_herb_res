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
  dx::Float64 = 1.0, dg::Float64 = 1.0, int_sd::Float64 = 1.4142, num_iter_pre::Int64 = 50,
  num_iter_post::Int64 = 50, lower_g:Float64 = -10.0, upper_g::Float64 = 10.0, offspring_sd:Float64 = 1.0, 
  germ_prob = 0.5, fec_0 = 10.0, fec_cost = 10.0, fec_max = 60.0, dd_fec = 0.1, 
  base_sur = 10.0, herb_effect = 20.0, g_prot = 10.0, seed_sur = 0.5, pro_exposed = 0.8,
  scale_pollen = 32.0, shape_pollen = 3.32, seed_pro_short = 0.48, seed_mean_short = 0.58,
  seeds_to_mean_short = 0.44, seed_mean_long = 1.65, seeds_to_mean_long = 0.39)
  
  # First the populations must be run for num_iter_pre to let the populations develope 
  seed_disp_mat_1D = zeros(convert(Int32, (length(int_rr_1) / dx) + 1), 
    convert(Int32, (length(int_rr_1) / dx) + 1))
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, seed_mean_dist_short, 
    pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long)
  
  pollen_disp_mat = zeros(convert(Int32, (length(int_rr_1) / dx) + 1), 
    convert(Int32, (length(int_rr_1) / dx) + 1))
  pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx, a = scale_pollen, c = shape_pollen)
   
  #set aside a chunck of memory for the landscapes for each genotype 
  RR_landscape = Array{Array{Float64, 2}, 1}(num_iter_pre + num_iter_post)
  Rr_landscape = Array{Array{Float64, 2}, 1}(num_iter_pre + num_iter_post)
  rr_landscape = Array{Array{Float64, 2}, 1}(num_iter_pre + num_iter_post)
  
  for t in 1:num_iter
    RR_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
    Rr_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
    rr_landscape[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
  end
  
  
  # TODO: the rest of the code
  
  
end