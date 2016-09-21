# simulation experiments for the TSR and NTSR joint evolution
# to test how TSR can invade into a population with NTSR

# Ideally I will be be able to specify which parameters to vary and
# the values to test them at. I could just pass every parameter 
# then I get length and run a for loop over it.

# need to run each source scenario three times, once in a empty landscape, one in an exposed population
# and one in a herbicide exposed population (can pre calculate all these so don't need to re-caclulate each time 

function transplant_exper(int_rr_1::Array{Float64, 1}; int_rr_2::Array{Float64, 1}, 
  int_Rr_1::Array{Float64, 1} = [1, 2], int_Rr_2::Array{Float64, 1} = [1, 2], 
  herb_app_source:Array{Float64, 1} = [1, 2], herb_app_pre_sink::Array{Float64, 1} = [1,2],
  herb_app_post_sink::Array{Float64, 1} = [1,2], dx::Float64 = 1.0, dg::Float64 = 1.0, 
  int_sd::Float64 = 1.4142, num_iter::Int64 = 50, lower_g:Float64 = -10.0, upper_g::Float64 = 10.0, 
  offspring_sd:Float64 = 1.0, germ_prob = 0.5, fec_0 = 10.0, fec_cost = 10.0, fec_max = 60.0, dd_fec = 0.1, 
  base_sur = 10.0, herb_effect = 20.0, g_prot = 10.0, seed_sur = 0.5, pro_exposed = 0.8,
  scale_pollen = 32.0, shape_pollen = 3.32, seed_pro_short = 0.48, seed_mean_short = 0.58,
  seeds_to_mean_short = 0.44, seed_mean_long = 1.65, seeds_to_mean_long = 0.39)
  
  # set up the seed and pollen dispersal kernels
  seed_disp_mat_1D = zeros(convert(Int32, (length(int_rr_1) / dx) + 1), 
    convert(Int32, (length(int_rr_1) / dx) + 1))
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, seed_mean_dist_short, 
    pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long)
  
  pollen_disp_mat = zeros(convert(Int32, (length(int_rr_1) / dx) + 1), 
    convert(Int32, (length(int_rr_1) / dx) + 1))
  pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx, a = scale_pollen, c = shape_pollen)
   
  #set aside a chunck of memory for the landscapes for each genotype in both the source and sink
  RR_landscape_source = Array{Array{Float64, 2}, 1}(num_iter)
  Rr_landscape_source = Array{Array{Float64, 2}, 1}(num_iter)
  rr_landscape_source = Array{Array{Float64, 2}, 1}(num_iter)
  
  RR_landscape_sink = Array{Array{Float64, 2}, 1}(num_iter)
  Rr_landscape_sink = Array{Array{Float64, 2}, 1}(num_iter)
  rr_landscape_sink = Array{Array{Float64, 2}, 1}(num_iter)
  
  for t in 1:num_iter
    RR_landscape_source[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
    Rr_landscape_source[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
    rr_landscape_source[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])  
    
    RR_landscape_sink[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
    Rr_landscape_sink[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
    rr_landscape_sink[t] = zeros(length(g_vals), size(seed_disp_mat_1D)[2])
  end
  
  # set an intial population in the source and sink and run for num_iter_pre
  
  # take a slice across g averaged over all locations in the last time step of source and add it to a single locaiton in sink
  
  # run sink for num_iter_post under herbicide application
  
  
end