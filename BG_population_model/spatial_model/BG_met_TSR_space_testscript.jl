# script to tes the functions produce expected results in simple cases 

using Base.Test
using Distributions

# set up some basic handelers to check fo rsuccesses, failuers and errors
cust_hand(r::Test.Success) = println("Success on $(r.expr)")
#throw error on failur so it shows the whole call tree
cust_hand(r::Test.Failure) = error("Failure on $(r.expr)")
cust_hand(r::Test.Error) = rethrow(r)

#pull in the functions
cd("/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/spatial_model")
include("BG_met_TSR_space_pop_process.jl")
include("BG_met_TSR_space_dispersal_functions.jl")

#set up all the variables I will need
int_pop_RR = (0, 0, [9, 10, 11]);
int_pop_Rr = (0.0001, 0, [9, 10, 11]);
int_pop_rr = (100.0, 0, [9, 10, 11]);
int_sd = 1.41;
num_iter = 10;
landscape_size = 30;
dx = 0.5;
lower_g = -10;
upper_g = 10;
germ_prob = 0.7;
offspring_sd = 1.0;
fec0 = 10.0;
fec_cost = 2;
fec_max = 100.0;
dd_fec = 0.004;
dg = 0.5;
base_sur = 10.0;
resist_G = ["RR", "Rr"];
herb_effect = 10.0;
g_prot = 1.0;
pro_exposed = 0.8;
seed_sur = 0.5

g_vals = lower_g : dg : upper_g


#test for seed bank up date
sb_now = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));
for i in [5, 6, 7]
  sb_now[:, i] = pdf(Normal(0, 1.41), g_vals) * 100; 
end
sb_next = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));
seedbank_update!(sb_next, sb_now, seed_sur);

#test for germination, both number of germinates and num left in seed bank 
sb_next2 = deepcopy(sb_next);
ag_plants = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));
new_plants!(ag_plants, sb_next2, germ_prob)

#test that g_mixing kernel cloumns intergrate to 1
g_mixing_kernel = zeros(length(g_vals), length(g_vals) ^ 2);
fill_g_mixing_kernel!(g_mixing_kernel, offspring_sd, g_vals)
g_mixing_index = generate_index_pairs(g_vals);

#test the dispersal matrix 
seed_disp_mat_1D = zeros(convert(Int32, (landscape_size / dx) + 1), 
  convert(Int32, (landscape_size / dx) + 1));
seed_disp_mat_builder_1D!(seed_disp_mat_1D, dx);
#test to check the edge rows and cols sum to 0.5
seed_disp_colsum = sum(seed_disp_mat_1D, 1);
seed_disp_rowsum = sum(seed_disp_mat_1D, 2);

pollen_disp_mat_1D = zeros(convert(Int32, (landscape_size / dx) + 1), 
  convert(Int32, (landscape_size / dx) + 1));
pollen_disp_mat_builder_1D!(pollen_disp_mat_1D, res = dx);
#test to check the edge rows and cols sum to 0.5
pollen_disp_colsum = sum(pollen_disp_mat_1D, 1);
pollen_disp_rowsum = sum(pollen_disp_mat_1D, 2);

#test total amount of pollen produced in the landscapoe is the expected amount after dispersal
#calculate the expected amount of polleen given the ag_plants and disperal kernel
#tot_pollen = sum((sum(ag_plants[:, [5, 6, 7]] , 1) * dg) .* sum(pollen_disp_mat_1D[:, [5, 6, 7]], 1))
tot_pollen_fixed = 89.24850144306271; #save so we have fixed reference point in case of future changes to functions

pollen_RR = zeros(length(g_vals), size(seed_disp_mat_1D)[1]);
pollen_Rr = zeros(length(g_vals), size(seed_disp_mat_1D)[1]);
pollen_rr = zeros(length(g_vals), size(seed_disp_mat_1D)[1]);
total_pollen = zeros(size(seed_disp_mat_1D)[1]);
## pollen at arrived each location for each g from all other locations  
pollen_RR[:, :] = ag_plants * pollen_disp_mat_1D * 0.1;
pollen_Rr[:, :] = ag_plants * pollen_disp_mat_1D * 0.2;
pollen_rr[:, :] = ag_plants * pollen_disp_mat_1D;
total_pollen[:] = sum(pollen_RR, 1) + sum(pollen_Rr, 1) + sum(pollen_rr, 1);
#normalise the pollen counts, double intergration across g and and all x as pollen comes from all x
pollen_RR[:, :] = pollen_RR ./ transpose(total_pollen * dg);
pollen_Rr[:, :] = pollen_Rr ./ transpose(total_pollen * dg);
pollen_rr[:, :] = pollen_rr ./ transpose(total_pollen * dg);

#test number of seeds produced at each location is expected
#se up new seed landscapes 
RR_newseed = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));
Rr_newseed = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));
rr_newseed = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));

g_effect_fec = exp(-(fec0 - g_vals * 0.0)); #first test assume that there is no fec cost and no effect of density

#work out how many seeds we expect from each G
expect_RR_seeds = sum(ag_plants * 0.1 .* ((fec_max / 3) ./ (1 + g_effect_fec)), 1) * dg;
expect_Rr_seeds = sum(ag_plants * 0.2 .* ((fec_max / 3) ./ (1 + g_effect_fec)), 1) * dg;
expect_rr_seeds = sum(ag_plants .* ((fec_max / 3) ./ (1 + g_effect_fec)), 1) * dg;

new_seeds_at_t!(RR_newseed, Rr_newseed, rr_newseed, ag_plants * 0.1, ag_plants * 0.2, ag_plants,
  pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index, g_effect_fec,
  fec_max, 0.0, dg)

#check the matrix-mult version produces the same answer
RR_newseed_mm = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));
Rr_newseed_mm = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));
rr_newseed_mm = zeros(length(g_vals), convert(Int32, (landscape_size / dx) + 1));
new_seeds_at_t_mm!(RR_newseed_mm, Rr_newseed_mm, rr_newseed_mm, ag_plants * 0.1, ag_plants * 0.2, ag_plants,
  pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index, g_effect_fec,
  fec_max, 0.0, dg)

# test the number of seeds dispersed to each location is the expected amount 
# the  number of expected seeds is the the number of seeds at each locaiton by the sum of the disp mat row for each of those locations
#num_exp = sum((sum(RR_newseed, 1) * dg) .* sum(seed_disp_mat_1D, 1))
tot_seed_fixed = 105.5967972459612;
tot_sb = sum(sb_next2) * dg;
sb_next3 = deepcopy(sb_next2);
sb_next3[:, :] = sb_next3  + (RR_newseed * seed_disp_mat_1D); 

#test the survival pre_calc works as expected 
sur_tup = survival_pre_calc(base_sur, convert(Array{Float64, 1}, g_vals), herb_effect, 
  g_prot, pro_exposed)
# sur_tup[1] == sur_tup[2] when (herb_effect - (g * g_prot)) == base_sur, which happen when g = 10 (g_vals[end])   
# also sur_tup[2][1] should be pretty close (but not exactly) pro_exposed
  
# expected number of survivours at each location  
herb_application = ones(size(ag_plants)[2]);
herb_application[[5, 7]] = 2;

ag_surs = deepcopy(ag_plants);
survival_at_t!(ag_surs, resist_G, "Rr", convert(Array{Int64, 1}, herb_application), sur_tup)
sur_resist = deepcopy(ag_surs[:, [5, 6, 7]]);


ag_surs = deepcopy(ag_plants);
survival_at_t!(ag_surs, resist_G, "rr", convert(Array{Int64, 1}, herb_application), sur_tup)
sur_sucep = deepcopy(ag_surs[:, [5, 6, 7]]);
#reference number so changes are picked up 
sur_herb_ref = 41.9993644298266; 

# test for the dist summary function that gets the mean, sd and total BG_population_model
test_dist = pdf(Normal(5, 1), g_vals) * 101
test_summary = dist_summary(test_dist, collect(g_vals), dg);

Test.with_handler(cust_hand) do
  #test the seedbank_update intergrates to 50 seeds
  a = sum(sb_next, 1) * dg
  @test isapprox(a[[5 6 7]], [50 50 50], atol = 0.000001)
  
  #test above and beloe groung germination
  ag = sum(ag_plants, 1) * dg
  @test isapprox(ag[[5 6 7]], [35 35 35], atol = 0.000001)
  sb = sum(sb_next2, 1) * dg
  @test isapprox(sb[[5 6 7]], [15 15 15], atol = 0.00001)
  
  #test the g mixing kernel intergrates to the expected values 
  @test isapprox(sum(g_mixing_kernel[:, [100, 101]], 1) * dg, [1.0 1.0], atol = 0.000001)
  edge_dist = pdf(Normal(lower_g, offspring_sd), g_vals)
  @test isapprox(sum(g_mixing_kernel[:, [1, end]], 1) * dg, [(sum(edge_dist) * dg) (sum(edge_dist) * dg)], atol = 0.00001)
  
  #test the disperal matrix is as expected 
  @test isapprox(seed_disp_colsum[[1 end]], [0.5 0.5], atol = 0.000001)
  @test seed_disp_colsum[convert(Int32, floor(length(seed_disp_colsum) / 2))] >= 0.999
  @test isapprox(seed_disp_rowsum[[1 end]], [0.5 0.5], atol = 0.000001)
  @test seed_disp_rowsum[convert(Int32, floor(length(seed_disp_colsum) / 2))] >= 0.999
  @test seed_disp_mat_1D[1, :] == transpose(seed_disp_mat_1D[:, 1])
  
  #pollen dispersal, more pollen lost off the edge of the landscape becasue the kernel is more spread
  @test isapprox(pollen_disp_colsum[[1 end]], [0.5 0.5], atol = 0.003)
  @test pollen_disp_colsum[convert(Int32, floor(length(pollen_disp_colsum) / 2))] >= 0.98
  @test isapprox(pollen_disp_rowsum[[1 end]], [0.5 0.5], atol = 0.003)
  @test pollen_disp_rowsum[convert(Int32, floor(length(pollen_disp_colsum) / 2))] >= 0.98
  @test pollen_disp_mat_1D[1, :] == transpose(pollen_disp_mat_1D[:, 1])
  @test isapprox(sum(total_pollen) * dx, tot_pollen_fixed * 1.3, atol = 0.000001) #1.3 casue 0.1 + 0.2 + 1 = 1.3, which are realtive sizes of ag_plants used to produce pollen 
  
  #should intergrate to 1 across g and G at each location
  @test isapprox((sum(pollen_RR, 1) + sum(pollen_Rr, 1) + sum(pollen_rr, 1)) * dg, 
    ones(1, size(total_pollen)[1]), atol = 0.000001)
  
  #test that the function produces the expected number of total seeds after mixing 
  @test isapprox((sum(rr_newseed, 1) + sum(RR_newseed, 1) + sum(Rr_newseed, 1))*dg, 
    sum(expect_rr_seeds, 1) + sum(expect_Rr_seeds, 1) + sum(expect_RR_seeds, 1), atol = 0.000001)
  @test isapprox(RR_newseed, RR_newseed_mm, rtol = 0.0000000000001)
  @test isapprox(Rr_newseed, Rr_newseed_mm, rtol = 0.0000000000001)
  @test isapprox(rr_newseed, rr_newseed_mm, rtol = 0.0000000000001)
  
  #test survival
  @test isapprox(sur_tup[1], sur_tup[2][end], atol = 0.0000000001)
  @test isapprox(sur_tup[2][1], (1 - pro_exposed), atol = 0.0001)
  @test isapprox(sum(sur_resist, 1)[1], sum(ag_plants, 1)[5] * (1 / (1 + exp(-base_sur))), atol = 0.000000001)
  @test isapprox(sum(sur_sucep, 1)[2], sum(sur_resist, 1)[1], atol = 0.0000000001)
  @test isapprox(sum(sur_sucep, 1)[1], sur_herb_ref, atol = 0.0000000001)
  
  #test seed dispersal 
  @test isapprox(sum(sb_next3) * dg * dx, (tot_sb + tot_seed_fixed) * dx, atol = 0.0000001)
  
  #test the dist summary function 
  @test isapprox(test_summary, [5.0, 1.0, 101.0], atol = 0.0001)

end


#BENCH MARKING 

#test the speed of matrix mult vs loop version of gene mixing

@time  new_seeds_at_t!(RR_newseed, Rr_newseed, rr_newseed, ag_plants * 0.1, ag_plants * 0.2, ag_plants,
  pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index, g_effect_fec,
  fec_max, 0.0, dg)

 @time  new_seeds_at_t!(RR_newseed, Rr_newseed, rr_newseed, ag_plants * 0.1, ag_plants * 0.2, ag_plants,
  pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index, g_effect_fec,
  fec_max, 0.0, dg)
  
 @time  new_seeds_at_t_mm!(RR_newseed, Rr_newseed, rr_newseed, ag_plants * 0.1, ag_plants * 0.2, ag_plants,
  pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index, g_effect_fec,
  fec_max, 0.0, dg)

 @time  new_seeds_at_t_mm!(RR_newseed, Rr_newseed, rr_newseed, ag_plants * 0.1, ag_plants * 0.2, ag_plants,
  pollen_RR, pollen_Rr, pollen_rr, g_mixing_kernel, g_mixing_index, g_effect_fec,
  fec_max, 0.0, dg)
 
 # Matrix-mult version about 1/3 faster and there are many fewer memory allocation 