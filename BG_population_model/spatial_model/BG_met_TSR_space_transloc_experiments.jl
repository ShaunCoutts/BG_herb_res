# simulation experiments for the TSR and NTSR joint evolution
# to test how TSR can invade into a population with NTSR

# Ideally I will be be able to specify which parameters to vary and
# the values to test them at. I could just pass every parameter 
# then I get length and run a for loop over it.

# need to run each source scenario three times, once in a empty landscape, one in an exposed population
# and one in a herbicide exposed population (can pre calculate all these so don't need to re-caclulate each time 
file_loc_func_p = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model" 
cd(file_loc_func_p)
include("BG_met_TSR_space_pop_process.jl")
include("BG_met_TSR_space_dispersal_functions.jl")
  

# single iteration takes a state of the population (i.e. one landscape with all TS genotypes) and progresses them
# 1 timestep foward, the first argument is the numerical next time step.
# the function takes the landscapes as references to arrays that can be modified in place (the xx_ls arrays), it all so takes a 
# set of holding arrays (xx_ab_pop, total_pollen, pollen_xx, xx_newseed, so these don't have to be built allocated every function call)

function one_step_foward!(next_t::Int64, RR_ls::Array{Float64, 3},  Rr_ls::Array{Float64, 3},
   rr_ls::Array{Float64, 3}, RR_ab_pop::Array{Float64, 2}, Rr_ab_pop::Array{Float64, 2}, 
   rr_ab_pop::Array{Float64, 2}, pollen_RR::Array{Float64, 2}, pollen_Rr::Array{Float64, 2},
   pollen_rr::Array{Float64, 2}, total_pollen::Array{Float64, 1}, RR_newseed::Array{Float64, 2}, 
   Rr_newseed::Arr{Float, 2}, rr_newseed::Array{Float64, 2}, dg::Float64, seed_sur::Float64, 
   germ_prob::Float64, sur_tup::Tuple{Float64, Array{Float64, 1}}, resist_G::Array{ASCIIString, 1},
   herb_application::Array{Int64, 1}, pollen_disp_mat::Array{Float64, 2}, fec_max::Float64, 
   dd_fec::Float64, g_mixing_kernel::Array{Float64, 2}, g_mixing_index::Array{Int64, 2}, 
   g_effect_fec::Array{Float64, 1})

  #move seeds to the next timestep, killing as we do so
  seedbank_update!(RR_ls[:, :, next_t], RR_ls[:, :, next_t - 1], seed_sur)
  seedbank_update!(Rr_ls[:, :, next_t], Rr_ls[:, :, next_t - 1], seed_sur)
  seedbank_update!(rr_ls[:, :, next_t], rr_ls[:, :, next_t - 1], seed_sur)
  
  #germination   
  new_plants!(RR_ab_pop, RR_ls[:, :, next_t], germ_prob)
  new_plants!(Rr_ab_pop, Rr_ls[:, :, next_t], germ_prob)
  new_plants!(rr_ab_pop, rr_ls[:, :, next_t], germ_prob)
  
  #above ground survival
  survival_at_t!(RR_ab_pop, resist_G, "RR", herb_application, sur_tup)
  survival_at_t!(Rr_ab_pop, resist_G, "Rr", herb_application, sur_tup)
  survival_at_t!(rr_ab_pop, resist_G, "rr", herb_application, sur_tup)
  
  ## pollen at arrived each location for each g from all other locations  
  pollen_RR[:, :] = RR_ab_pop * pollen_disp_mat
  pollen_Rr[:, :] = Rr_ab_pop * pollen_disp_mat
  pollen_rr[:, :] = rr_ab_pop * pollen_disp_mat
  total_pollen[:] = sum(pollen_RR, 1) + sum(pollen_Rr, 1) + sum(pollen_rr, 1)
  #normalise the pollen counts, double intergration across g, intergration across x already done in building dispersal matrix
  pollen_RR[:, :] = pollen_RR ./ transpose(total_pollen * dg)
  pollen_Rr[:, :] = pollen_Rr ./ transpose(total_pollen * dg)
  pollen_rr[:, :] = pollen_rr ./ transpose(total_pollen * dg)
  
  #create new seeds
  new_seeds_at_t_mm!(RR_newseed, Rr_newseed, rr_newseed, RR_ab_pop, Rr_ab_pop,
    rr_ab_pop, pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index,   
    g_effect_fec, fec_max, dd_fec, dg)
  
  # disperse the seeds
  RR_ls[:, :, next_t] = RR_ls[:, :, next_t]  + (RR_newseed * seed_disp_mat_1D)
  Rr_ls[:, :, next_t] = Rr_ls[:, :, next_t]  + (Rr_newseed * seed_disp_mat_1D)
  rr_ls[:, :, next_t] = rr_landscape[:, :, next_t]  + (rr_newseed * seed_disp_mat_1D)

  return nothing
  
end

# injects num_seeds into landscape_xx at locations locs, at time t_inject,
# with a given mean_g and sd_g.
function hot_seed_injection!(RR_ls_empty::Array{Float64, 3}, Rr_ls_empty::Array{Float64, 3}, 
  rr_ls_empty::Array{Float64, 3}, RR_ls_naive::Array{Float64, 3}, Rr_ls_naive::Array{Float64, 3}, 
  rr_ls_naive::Array{Float64, 3}, RR_ls_expos::Array{Float64, 3}, Rr_ls_expos::Array{Float64, 3}, 
  rr_ls_expos::Array{Float64, 3}, num_RR::Float64, num_Rr::Float64, num_rr::Float64, 
  locs::Array{Float64, 1}, t_inject::Int64, inject_mean_g::Float64, inject_sd_g::Float64, 
  g_vals::Array{Float64, 1})
  
  RR_ls_empty[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_RR
  Rr_ls_empty[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_Rr
  rr_ls_empty[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_rr
  
  RR_ls_naive[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_RR
  Rr_ls_naive[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_Rr
  rr_ls_naive[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_rr
  
  RR_ls_expos[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_RR
  Rr_ls_expos[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_Rr
  rr_ls_expos[:, locs, t_inject] += pdf(Normal(inject_mean_g, inject_sd_g), g_vals) * num_rr
  
  return nothing

end

# Takes the snapshots of the population at each timestep based on the three G landscapes for a single scenario 
# returns a 11 x time_steps 2D array, first 4 rows give total pops (3 G types + total), next 4 mean_g (3 G types + total), 
# next 3 %occ (3 G types) 
function pop_snapshots(RR_ls::Array{Float64, 3}, Rr_ls::Array{Float64, 3},
  rr_ls::Array{Float64, 3}, g_vals::Array{Float64, 1}, dg::Float64, dx::Float64)
  
  time_steps = size(RR_ls)[3]
  output = zeros(11, time_steps) # 11 because #G (3), mean_g for each G (3), mean_g overall, total_pop, #occ each G (3) (3 + 3 + 3 + 1 + 1 = 11)  
  
  # total populations for each G
  output[1, :] = sum(RR_ls, [1, 2]) * dg * dx
  output[2, :] = sum(Rr_ls, [1, 2]) * dg * dx
  output[3, :] = sum(rr_ls, [1, 2]) * dg * dx
  # total population
  output[4, :] = sum(output[1:3, :], 1)
 
  # mean_g
  output[5, :] = (sum((sum(RR_ls, 2) * dx) .* g_vals, 1) * dg) ./ output[1, :]
  output[6, :] = (sum((sum(Rr_ls, 2) * dx) .* g_vals, 1) * dg) ./ output[2, :]
  output[7, :] = (sum((sum(rr_ls, 2) * dx) .* g_vals, 1) * dg) ./ output[3, :]
  # population mean_g 
  output[8, :] = (sum(((sum(RR_ls, 2) + sum(Rr_ls, 2) + sum(rr_ls, 2)) * dx) .* g_vals, 1) * dg) ./ output[4, :]
  
  # % landscape occupied 
  output[9, :] = sum(sum(RR_ls, 1) * dg .> 1.0, 2) / size(RR_ls)[2]
  output[10, :] = sum(sum(Rr_ls, 1) * dg .> 1.0, 2) / size(Rr_ls)[2]
  output[11, :] = sum(sum(rr_ls, 1) * dg .> 1.0, 2) / size(rr_ls)[2]
  
  return output
  
end

# run an entire sceanrio to produce a set of data for a plot of three measures of populaiton 
# resistance and extent over time (enough data for 9 plots, 3 measures x 3 sink pop types)
function run_scene_trans(g_vals::Array{Float64, 1}, x_dim::Int64, dg::Float64, dx::Float64, 
  num_iter::Int64, burnin::Int64, num_inject::Float64, pro_R_inject::Float64, inject_mean_g::Float64, 
  inject_sd_g::Float64, inject_locs::Array{Float64, 1}, int_rr::Float64, int_mean_g::Float64,
  int_sd_g::Float64,
  seed_sur::Float64, germ_prob::Float64, resist_G::Array{ASCIIString, 1}, fec_max::Float64, 
  dd_fec::Float64, fec0::Float64, fec_cost::Float64, base_sur::Float64, herb_effect::Float64, 
  g_prot::Float64, pro_exposed::Float64,
  
  seed_pro_short::Float64, seed_mean_dist_shor::Float64, pro_seeds_to_mean_short::Float64, 
  seed_mean_dist_long::Float64, pro_seeds_to_mean_long::Float64, scale_pollen::Float64, shape_pollen::Float64,
  offspring_sd::Float64, 
  )

  # set of temporary holdoing matricies to set aside some memory, so these don't have to be rebuilt at each iteration
  # a set of matrices to hold the above ground populations
  RR_ab_pop = zeros(length(g_vals), x_dim)
  Rr_ab_pop = zeros(length(g_vals), x_dim)
  rr_ab_pop = zeros(length(g_vals), x_dim)

  # a set of matrices to hold the total amount of pollen that arrives are each location for each metabolic 
  # resitance score for each genotype
  pollen_RR = zeros(length(g_vals), x_dim)
  pollen_Rr = zeros(length(g_vals), x_dim)
  pollen_rr = zeros(length(g_vals), x_dim)
  total_pollen = zeros(x_dim)

  #set of matrices to hold the new seeds produced at each location pre dispersal 
  RR_newseed = zeros(length(g_vals), x_dim)
  Rr_newseed = zeros(length(g_vals), x_dim)
  rr_newseed = zeros(length(g_vals), x_dim)

  # set up the seed and pollen dispersal kernels 
  seed_disp_mat_1D = zeros(x_dim, x_dim)
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx, seed_pro_short, seed_mean_dist_short, 
    pro_seeds_to_mean_short, seed_mean_dist_long, pro_seeds_to_mean_long)

  pollen_disp_mat = zeros(x_dim, x_dim)
  pollen_disp_mat_builder_1D!(pollen_disp_mat, res = dx, a = scale_pollen, c = shape_pollen)

  # build the mixing kernel for metabolic resistance score every row is a offspring score
  # every coloum is a g_m x g_p combination, so every coloumn is a normal dist with
  # a mean of g_m*g_p and sd of offspring sd
  g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2)
  fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
  g_mixing_index = generate_index_pairs(g_vals)
  # give the effect of herb as a function of g, make it symetrical stabilising function, centered on 0
  g_effect_fec = exp(-(fec0 - abs(g_vals) * fec_cost))

  # set up the survival vectors 
  sur_tup = survival_pre_calc(base_sur, g_vals, herb_effect, g_prot, pro_exposed)

  # set up the three sink landscapes, empty, naive and exposed      
  # set aside a chunck of memory for the landscapes for each genotype or structure [g_vals, x_dim, timesteps] 
  RR_empty = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_empty = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_empty = zeros(length(g_vals), x_dim, num_iter + burnin)

  RR_naive = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_naive = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_naive = zeros(length(g_vals), x_dim, num_iter + burnin)

  RR_expos = zeros(length(g_vals), x_dim, num_iter + burnin)
  Rr_expos = zeros(length(g_vals), x_dim, num_iter + burnin)
  rr_expos = zeros(length(g_vals), x_dim, num_iter + burnin)

  # intitilase the population on the full landscapes in the first time slice
  for x in 1:x_dim
    rr_naive[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_rr
    rr_expos[:, x, 1] = pdf(Normal(int_mean_g, int_sd_g), g_vals) * int_rr
  end 
  
  #herb application for navie and exposed populations
  herb_app_naive = ones(x_dim)
  herb_app_expos = deepcopy(herb_app_naive)
  herb_app_expos = += 1 #adds one to the locations where herbicide is going to be applied 
  
  # run the burnin period on the exposed and niave pops
  for t in 2:burnin
  
    # step through the naive population, letting it develope
    one_step_foward!(t, RR_naive,  Rr_naive, rr_naive, RR_ab_pop, Rr_ab_pop, 
      rr_ab_pop, pollen_RR, pollen_Rr, pollen_rr, total_pollen, RR_newseed, 
      Rr_newseed, rr_newseed, dg, seed_sur, germ_prob, sur_tup, resist_G,
      herb_app_naive, pollen_disp_mat, fec_max, dd_fec, g_mixing_kernel, 
      g_mixing_index, g_effect_fec)
      
    # step through the exposed population letting it develope
    one_step_foward!(t, RR_expos,  Rr_expos, rr_expos, RR_ab_pop, Rr_ab_pop, 
      rr_ab_pop, pollen_RR, pollen_Rr, pollen_rr, total_pollen, RR_newseed, 
      Rr_newseed, rr_newseed, dg, seed_sur, germ_prob, sur_tup, resist_G,
      herb_app_expos, pollen_disp_mat, fec_max, dd_fec, g_mixing_kernel, 
      g_mixing_index, g_effect_fec)
      
  end
  
  # introduce seeds from the population with different TSR and quantitative resistance
  num_RR_inject = num_inject * pro_R_inject
  num_Rr_inject = 0.0
  num_rr_inject = 1 - num_RR_inject
  
  hot_seed_injection!(RR_ls_empty, Rr_ls_empty, rr_ls_empty, RR_ls_naive, Rr_ls_naive, 
    rr_ls_naive, RR_ls_expos, Rr_ls_expos, rr_ls_expos, num_RR_inject, num_Rr_inject, 
    num_rr_inject, inject_locs, burnin + 1, inject_mean_g, inject_sd_g, g_vals)
  
  # run the model for num_iter under herbicide application  
  for t in (burnin + 1):(num_iter + burnin)
  
    # step through the naive population
    one_step_foward!(t, RR_naive,  Rr_naive, rr_naive, RR_ab_pop, Rr_ab_pop, 
      rr_ab_pop, pollen_RR, pollen_Rr, pollen_rr, total_pollen, RR_newseed, 
      Rr_newseed, rr_newseed, dg, seed_sur, germ_prob, sur_tup, resist_G,
      herb_app_expos, pollen_disp_mat, fec_max, dd_fec, g_mixing_kernel, 
      g_mixing_index, g_effect_fec)
      
    # step through the exposed population
    one_step_foward!(t, RR_expos,  Rr_expos, rr_expos, RR_ab_pop, Rr_ab_pop, 
      rr_ab_pop, pollen_RR, pollen_Rr, pollen_rr, total_pollen, RR_newseed, 
      Rr_newseed, rr_newseed, dg, seed_sur, germ_prob, sur_tup, resist_G,
      herb_app_expos, pollen_disp_mat, fec_max, dd_fec, g_mixing_kernel, 
      g_mixing_index, g_effect_fec)
     
    # step through the empty population
    one_step_foward!(t, RR_empty,  Rr_empty, rr_empty, RR_ab_pop, Rr_ab_pop, 
      rr_ab_pop, pollen_RR, pollen_Rr, pollen_rr, total_pollen, RR_newseed, 
      Rr_newseed, rr_newseed, dg, seed_sur, germ_prob, sur_tup, resist_G,
      herb_app_expos, pollen_disp_mat, fec_max, dd_fec, g_mixing_kernel, 
      g_mixing_index, g_effect_fec)
      
  end
  
pop_snapshots(RR_ls::Array{Float64, 3}, Rr_ls::Array{Float64, 3},
  rr_ls::Array{Float64, 3}, g_vals::Array{Float64, 1}, dg::Float64, dx::Float64)

  output_empty = pop_snapshots(RR_empty, Rr_empty, rr_empty, g_vals, dg, dx)
  output_naive = pop_snapshots(RR_naive, Rr_naive, rr_naive, g_vals, dg, dx)
  output_expos = pop_snapshots(RR_expos, Rr_expos, rr_expos, g_vals, dg, dx)
  
  return (output_empty, output_naive, output_expos)
  
end

###################################################################################################
# script to run the translocation experiments
# parameter values
upper_g = 10;
lower_g = -10;
dg = 0.5;
int_g_mean = 0.0;
int_g_sd = 1.4142;

x_dim = 500; # number of spatial evaluation points, actual landscape size is x_dim * dx
dx = 1.0;

int_num_RR = 0;
int_num_Rr = 0;
int_num_rr = 10; # number of intial seeds at each location for each genoptype, assume only TS susceptible
burnin = 20;
num_iter = 50;

seed_pro_short = 0.4; 
seed_mean_dist_short = 0.5; 
pro_seeds_to_mean_short = 0.4; 
seed_mean_dist_long = 1.5; 
pro_seeds_to_mean_long = 0.4;
scale_pollen = 32.0;
shape_pollen = 3.23; 
offspring_sd = 1.0;
fec0 = 10.0;
fec_cost = 1.0;
base_sur = 10.0; 
herb_effect = 20.0; 
g_prot = 1.0; 
pro_exposed = 0.8;


# set up the evaluation points for quantitative resistance
g_vals = collect(lower_g : dg : upper_g);   















