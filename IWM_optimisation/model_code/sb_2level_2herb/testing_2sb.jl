# testing script to make sure the population model gives sensible results 
using Base.Test

#pull in the functions
cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb/")
include("pop_process_2sb_2herb.jl")
include("managment_functions.jl")


cov_mat = [1.0 0.0;
	   0.0 1.0]
dg = 2.0;
g_vals = collect(-12 : dg : 12);
len_g = size(g_vals)[1];
g1_vals = repeat(g_vals, inner = len_g);
g2_vals = repeat(g_vals, outer = len_g);

# define parameters 
fr = 0.5;
f0 = 4.0;
new_seeds = zeros(size(g1_vals)[1]);

pro_exposed = 0.8;
sur_base = 10.0;
effect_herb1 = 15.0; 
effect_herb2 = 15.0;
prot_g1_herb1 = 2.0; 
prot_g1_herb2 = 1.0;
prot_g2_herb1 = 1.0;
prot_g2_herb2 = 2.0;
sur_crop_alt = 0.8;
time_horizion = 20;
inv_frac = 0.8;
germ_prob = 0.4;
seed_sur = 0.6;
fec_max = 60.0; 
fec_dd = 0.004;

  
# define precomputed objects
crop_sur_tup = (1.0, sur_crop_alt, 0.0);

herb_sur_tup = survial_herb_setup(g1_vals, g2_vals, pro_exposed, sur_base, 
  effect_herb1, effect_herb2, prot_g1_herb1, prot_g1_herb2, prot_g2_herb1, prot_g2_herb2);

plow_subact = (false, true);

mix_keys = make_index_keys(len_g, len_g);

off_tem = offspring_dist_setup(g1_vals, g2_vals, cov_mat, mix_keys);

fec_cost = fec_cost_maker(fr, f0, g1_vals, g2_vals);

# define holding arrays
seedbank_1 = zeros(time_horizion, size(g1_vals)[1]);
seedbank_2 = zeros(time_horizion, size(g1_vals)[1]);
mat_pop = zeros(time_horizion, size(g1_vals)[1]);
pat_pop = zeros(size(g1_vals)[1]);
new_seeds = zeros(size(g1_vals)[1]);
par_mix = zeros(size(mix_keys[1]));

action_space = make_action_space();
act_seq = act_seq_herb0(action_space, time_horizion)

#make the plow array
plow_tup = plow_subact[act_seq[:, ACT_PLOW]];
plow_seq = collect(plow_tup);

# create an intial population 
int_cov = [1.414 0.0; 0.0 1.414];
int_dist = MvNormal(int_cov);
int_sb1 = pdf(int_dist, transpose(hcat(g1_vals, g2_vals))); 
int_sb2 = pdf(int_dist, transpose(hcat(g1_vals, g2_vals)));

# test of germination seed mort function
t = 1;
sb1_test = deepcopy(seedbank_1);
sb1_test[t, :] = int_sb1 * 10.0;
tot_seeds_t = sum(sb1_test[t, :]) * dg * dg;
sb2_test = deepcopy(seedbank_2);
sb2_test = sb2_test * 0.0;

germ_seed_sur!(sb1_test, sb2_test, mat_pop, t, germ_prob, seed_sur);

obs_sm_test = sum(sb1_test[t, :]) * dg * dg;
exp_sm_test = tot_seeds_t * seed_sur * (1 - germ_prob);

obs_sm_test_ab = sum(mat_pop[t, :]) * dg * dg;
exp_sm_test_ab = tot_seeds_t * seed_sur * germ_prob;

# test seed invertion 
t = 1;
sb1_test = deepcopy(seedbank_1);
sb1_test[t, :] = int_sb1 * 10.0;
tot_seeds_sb1 = sum(sb1_test[t, :]) * dg * dg;
sb2_test = deepcopy(seedbank_2);
sb2_test[t, :] = int_sb2 * 20.0;
tot_seeds_sb2 = sum(sb2_test[t, :]) * dg * dg;

plow_inversion!(sb1_test, sb2_test, t, true, inv_frac);

obs_plow_test_sb1 = sum(sb1_test[t, :]) * dg * dg;
obs_plow_test_sb2 = sum(sb2_test[t, :]) * dg * dg;

exp_plow_sb1 = tot_seeds_sb1 - (tot_seeds_sb1 * inv_frac) + (tot_seeds_sb2 * inv_frac); 
exp_plow_sb2 = tot_seeds_sb2 - (tot_seeds_sb2 * inv_frac) + (tot_seeds_sb1 * inv_frac); 

exp_plow_tot = tot_seeds_sb1 + tot_seeds_sb2;

# test survival function
t = 1;
mat_pop[t, :] = deepcopy(int_sb1) * 10.0;
tot_pre_sur  = sum(mat_pop[t, :]) * dg * dg;

survivours_herb0  = (mat_pop[t, :] .* herb_sur_tup[HERB0]) * crop_sur_tup[CROP_WW];
survivours_herb1  = (mat_pop[t, :] .* herb_sur_tup[HERB1]) * crop_sur_tup[CROP_WW];
survivours_herb2  = (mat_pop[t, :] .* herb_sur_tup[HERB2]) * crop_sur_tup[CROP_WW];
survivours_herb12  = (mat_pop[t, :] .* herb_sur_tup[HERB12]) * crop_sur_tup[CROP_WW];
# fix the num sur expected as a reference 
# exp_num_sur_herb1 = sum(survivours_herb1) * dg * dg;
exp_num_sur_herb1 = 2.4988759489825263;
exp_num_sur_herb12 = 2.155935161260673; 

# test seed production
t = 1;
mat_pop[t, :] = deepcopy(int_sb1) * 10.0;
tot_ab  = sum(mat_pop[t, :]) * dg * dg;
new_seeds = zeros(size(g1_vals)[1]);

fec_cost = fec_cost_maker(fr, f0, g1_vals, g2_vals);
make_new_seeds!(mat_pop, t, new_seeds, fec_max, tot_ab, fec_dd, fec_cost);

# exp_num_seeds = sum(new_seeds) * dg * dg;
exp_num_seeds = 554.1168280366774;
obs_num_newseed = sum(new_seeds) * dg * dg;

# decrease fec_cost to 0
fec_cost = fec_cost_maker(0.0, 10.0, g1_vals, g2_vals);
make_new_seeds!(mat_pop, t, new_seeds, fec_max, tot_ab, 0.0004, fec_cost);
seed_num_no_cost = sum(new_seeds) * dg * dg

# test seed mixing that there is the same number of seeds before and after mixing
# set it up so seeds are only produced in the center to stop eviction
seed_mat = zeros(size(g1_vals)[1]);
seed_mat[85] = 10.0; # put 10 seeds in the center of the (g1, g2 space) 
pat_mat = zeros(size(g1_vals)[1]);
pat_mat[85] = 1.0;
pat_mat = pat_mat / (sum(pat_mat) * dg * dg);


cs = seed_cross(pat_mat, seed_mat, dg, par_mix, mix_keys, off_tem);


@testset "all_tests" begin
  @testset "pop_model_tests" begin
  
    @testset "seed mort and germ" begin 
      @test isapprox(obs_sm_test, exp_sm_test, atol = 0.001);
      @test isapprox(obs_sm_test_ab, exp_sm_test_ab, atol = 0.001);
    end
  
    @testset "plowing" begin 
      @test isapprox(obs_plow_test_sb1, exp_plow_sb1, atol = 0.001);
      @test isapprox(obs_plow_test_sb2, exp_plow_sb2, atol = 0.001);
      @test isapprox(obs_plow_test_sb1 + obs_plow_test_sb2, exp_plow_tot, atol = 0.001);
    end
    
    @testset "survival" begin
      @test isapprox(sum(survivours_herb0) * dg * dg, tot_pre_sur, atol = 0.1);
      @test isapprox(sum(survivours_herb1) * dg * dg, exp_num_sur_herb1, atol = 0.00001);
      @test isapprox(sum(survivours_herb2) * dg * dg, exp_num_sur_herb1, atol = 0.00001);
      @test isapprox(sum(survivours_herb12) * dg * dg, exp_num_sur_herb12, atol = 0.00001);
    end
    
    @testset "seed production and mixing" begin
      @test isapprox(obs_num_newseed, exp_num_seeds, atol = 0.00001);
      @test isapprox(seed_num_no_cost, tot_ab * fec_max, atol = 5.0);
      @test isapprox(sum(cs) * dg * dg, sum(seed_mat) * dg * dg, atol = 3.0) # note there is about a 3% inaccuracy when dg = 2.0
    end
 
  end
  
end

################################################################################################################################3333
# I also want to make sure the selection is working as expected, make an animated heatmap to show how the population changes 

using Plots
pyplot()

cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb/")
include("pop_process_2sb_2herb.jl")
include("managment_functions.jl")


cov_mat = [1.0 0.0;
	   0.0 1.0]
dg = 2.0;
g_vals = collect(-12 : dg : 12);
len_g = size(g_vals)[1];
g1_vals = repeat(g_vals, inner = len_g);
g2_vals = repeat(g_vals, outer = len_g);

# define parameters 
fr = 0.5;
f0 = 4.0;
new_seeds = zeros(size(g1_vals)[1]);

pro_exposed = 0.8;
sur_base = 10.0;
effect_herb1 = 15.0; 
effect_herb2 = 15.0;
prot_g1_herb1 = 2.0; 
prot_g1_herb2 = 1.0;
prot_g2_herb1 = 1.0;
prot_g2_herb2 = 2.0;
sur_crop_alt = 0.8;
time_horizion = 25;
inv_frac = 0.8;
germ_prob = 0.4;
seed_sur = 0.6;
fec_max = 60.0; 
fec_dd = 0.004;

  
# define holding arrays
sb1_herb0_ncr = zeros(time_horizion, size(g1_vals)[1]);
sb1_herb1_ncr = zeros(time_horizion, size(g1_vals)[1]);
sb1_herb2_ncr = zeros(time_horizion, size(g1_vals)[1]);
sb1_herb12_ncr = zeros(time_horizion, size(g1_vals)[1]);
sb1_herb0_cr = zeros(time_horizion, size(g1_vals)[1]);
sb1_herb1_cr = zeros(time_horizion, size(g1_vals)[1]);
sb1_herb2_cr = zeros(time_horizion, size(g1_vals)[1]);
sb1_herb12_cr = zeros(time_horizion, size(g1_vals)[1]);

seedbank_2 = zeros(time_horizion, size(g1_vals)[1]);
mat_pop = zeros(time_horizion, size(g1_vals)[1]);
pat_pop = zeros(size(g1_vals)[1]);
new_seeds = zeros(size(g1_vals)[1]);

action_space = make_action_space();
act_seq_h0 = act_seq_herb0(action_space, time_horizion);
act_seq_h1 = act_seq_herb1(action_space, time_horizion);
act_seq_h2 = act_seq_herb2(action_space, time_horizion);
act_seq_h12 = act_seq_herb12(action_space, time_horizion);

#make the plow array
plow_subact = (false, true);
plow_tup = plow_subact[act_seq_h0[:, ACT_PLOW]];
plow_seq = collect(plow_tup);

# create an intial population 
int_cov = [1.414 0.0; 0.0 1.414];
int_dist = MvNormal(int_cov);
int_sb1 = pdf(int_dist, transpose(hcat(g1_vals, g2_vals))); 
int_sb2 = pdf(int_dist, transpose(hcat(g1_vals, g2_vals)));

# define precomputed objects
crop_sur_tup = (1.0, sur_crop_alt, 0.0);

herb_sur_tup_cr = survial_herb_setup(g1_vals, g2_vals, pro_exposed, sur_base, 
  effect_herb1, effect_herb2, prot_g1_herb1, prot_g1_herb2, prot_g2_herb1, prot_g2_herb2);

herb_sur_tup_ncr = survial_herb_setup(g1_vals, g2_vals, pro_exposed, sur_base, 
  effect_herb1, effect_herb2, prot_g1_herb1, 0.0, 0.0, prot_g2_herb2);
  

mix_keys = make_index_keys(len_g, len_g);

off_tem = offspring_dist_setup(g1_vals, g2_vals, cov_mat, mix_keys);

fec_cost = fec_cost_maker(fr, f0, g1_vals, g2_vals);

par_mix = zeros(size(mix_keys[1]));

#do the runs 
one_run!(sb1_herb0_ncr, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup_ncr, act_seq_h0[:, ACT_CROP], plow_seq, act_seq_h0[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg);
  
one_run!(sb1_herb1_ncr, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup_ncr, act_seq_h1[:, ACT_CROP], plow_seq, act_seq_h1[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg);
    
one_run!(sb1_herb2_ncr, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup_ncr, act_seq_h2[:, ACT_CROP], plow_seq, act_seq_h2[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg);
    
one_run!(sb1_herb12_ncr, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup_ncr, act_seq_h12[:, ACT_CROP], plow_seq, act_seq_h12[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg);
  
one_run!(sb1_herb0_cr, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup_cr, act_seq_h0[:, ACT_CROP], plow_seq, act_seq_h0[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg);
  
one_run!(sb1_herb1_cr, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup_cr, act_seq_h1[:, ACT_CROP], plow_seq, act_seq_h1[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg);
    
one_run!(sb1_herb2_cr, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup_cr, act_seq_h2[:, ACT_CROP], plow_seq, act_seq_h2[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg);
    
one_run!(sb1_herb12_cr, seedbank_2, mat_pop, pat_pop, new_seeds, par_mix, mix_keys, 
  off_tem, fec_cost, crop_sur_tup, herb_sur_tup_cr, act_seq_h12[:, ACT_CROP], plow_seq, act_seq_h12[:, ACT_HERB], 
  time_horizion, int_sb1, int_sb2, inv_frac, germ_prob, seed_sur, fec_max, fec_dd, dg);
 
 
 
 # set up the fonts fo all the labels first 
ax_font = Plots.Font("FreeSans", 12, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0));
title_font = Plots.Font("FreeSans", 14, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0));
leg_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0));
tic_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0));

grid_size = size(g_vals)[1];
adjust_scale = 1.1;
grey_pal = colormap("Grays", 100);
# make the animation object
anim = @animate for i = 1:time_horizion
  
  print(i)
  
  layout_arr = @layout grid(2, 4);

  plt = plot(layout = layout_arr, grid = true, background_color_outside = grey_pal[10], 
    border = false, foreground_color_grid = grey_pal[1],    
    title = ["herb 0: no change" "herb 1: move left" "herb 2: move up" "herb 1-2: move up 45 deg" "" "" "" ""],
    titleloc = :left, titlefont = title_font, guidefont = ax_font, tickfont = tic_font, 
    size = (600 * adjust_scale, 550 * adjust_scale),  
    xlabel = ["" "" "" "" "g1" "g1" "g1" "g1"], ylabel = ["no cross-res\ng2" "" "" "" "cross-res\ng2" "" "" ""]);

  plot_mat = reshape(sb1_herb0_ncr[i, :], (grid_size, grid_size));
  heatmap!(plt, plot_mat[5:end, 5:end], subplot = 1);
  
  plot_mat = reshape(sb1_herb1_ncr[i, :], (grid_size, grid_size));
  heatmap!(plt, plot_mat[5:end, 5:end], subplot = 2);
  
  plot_mat = reshape(sb1_herb2_ncr[i, :], (grid_size, grid_size));
  heatmap!(plt, plot_mat[5:end, 5:end], subplot = 3);

  plot_mat = reshape(sb1_herb12_ncr[i, :], (grid_size, grid_size));
  heatmap!(plt, plot_mat[5:end, 5:end], subplot = 4);
  
  plot_mat = reshape(sb1_herb0_cr[i, :], (grid_size, grid_size));
  heatmap!(plt, plot_mat[5:end, 5:end], subplot = 5);
  
  plot_mat = reshape(sb1_herb1_cr[i, :], (grid_size, grid_size));
  heatmap!(plt, plot_mat[5:end, 5:end], subplot = 6);
  
  plot_mat = reshape(sb1_herb2_cr[i, :], (grid_size, grid_size));
  heatmap!(plt, plot_mat[5:end, 5:end], subplot = 7);

  plot_mat = reshape(sb1_herb12_cr[i, :], (grid_size, grid_size));
  heatmap!(plt, plot_mat[5:end, 5:end], subplot = 8);
 
  plt
 
end 
 
# make the gif 
cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/");
gif(anim, "herb_effect_sb1.gif", fps = 2)
 
 
# above does not work very well, I will try a different approach 

adjust_scale = 2.0;

anim = @animate for i = 1:time_horizion
  
  print(i)
  
  layout_arr = @layout grid(2, 4);

  plt = plot(layout = layout_arr, grid = true, background_color_outside = grey_pal[10], 
    border = false, foreground_color_grid = grey_pal[1],    
    title = ["herb 0" "herb 1" "herb 2" "herb 1-2" "" "" "" ""],
    titleloc = :left, titlefont = title_font, guidefont = ax_font, tickfont = tic_font, 
    size = (600 * adjust_scale, 550 * adjust_scale),  
    xlabel = ["" "" "" "" "g1" "g1" "g1" "g1"], ylabel = ["no cross-res\ng2" "" "" "" "cross-res\ng2" "" "" ""]);

  plot_mat = reshape(sb1_herb0_ncr[i, :], (grid_size, grid_size));
  dist_g1 = vcat(sum(plot_mat, 1)...);
  dist_g2 = sum(plot_mat, 2);
  plot!(plt, [dist_g1 dist_g2], subplot = 1);
  
  plot_mat = reshape(sb1_herb1_ncr[i, :], (grid_size, grid_size));
  dist_g1 = vcat(sum(plot_mat, 1)...);
  dist_g2 = sum(plot_mat, 2);
  plot!(plt, [dist_g1 dist_g2], subplot = 2);
  
  plot_mat = reshape(sb1_herb2_ncr[i, :], (grid_size, grid_size));
  dist_g1 = vcat(sum(plot_mat, 1)...);
  dist_g2 = sum(plot_mat, 2);
  plot!(plt, [dist_g1 dist_g2], subplot = 3);

  plot_mat = reshape(sb1_herb12_ncr[i, :], (grid_size, grid_size));
  dist_g1 = vcat(sum(plot_mat, 1)...);
  dist_g2 = sum(plot_mat, 2);
  plot!(plt, [dist_g1 dist_g2], subplot = 4);
  
  plot_mat = reshape(sb1_herb0_cr[i, :], (grid_size, grid_size));
  dist_g1 = vcat(sum(plot_mat, 1)...);
  dist_g2 = sum(plot_mat, 2);
  plot!(plt, [dist_g1 dist_g2], subplot = 5);
  
  plot_mat = reshape(sb1_herb1_cr[i, :], (grid_size, grid_size));
  dist_g1 = vcat(sum(plot_mat, 1)...);
  dist_g2 = sum(plot_mat, 2);
  plot!(plt, [dist_g1 dist_g2], subplot = 6);
  
  plot_mat = reshape(sb1_herb2_cr[i, :], (grid_size, grid_size));
  dist_g1 = vcat(sum(plot_mat, 1)...);
  dist_g2 = sum(plot_mat, 2);
  plot!(plt, [dist_g1 dist_g2], subplot = 7);

  plot_mat = reshape(sb1_herb12_cr[i, :], (grid_size, grid_size));
  dist_g1 = vcat(sum(plot_mat, 1)...);
  dist_g2 = sum(plot_mat, 2);
  plot!(plt, [dist_g1 dist_g2], subplot = 8);
 
  plt
 
end 
 
cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/");
gif(anim, "herb_effect_sb1_lines.gif", fps = 2)

# all looks good