# script to tes the functions produce expected results in simple cases

using Base.Test

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
landscape_size = 20;
dx = 0.5;
lower_g = -10;
upper_g = 10;
germ_prob = 0.7;
offspring_sd = 1.0;
fec0 = 10;
fec_cost = 2;
fec_max = 100.0;
dd_fec = 0.004;
dg = 0.5;
base_sur = 10.0;
resist_G = ["RR", "Rr"];
herb_application = zeros(10);
herb_effect = 20.0;
g_prot = 1.0;
pro_exposed = 0.8;

#test for seed bank up date


#test for germination, both number of germinates and num left in seed bank 

#test that g_mixing kernel cloumns intergrate to 1

#test number of seeds produced at each location is expected

#test the dispersal matrix 
  seed_disp_mat_1D = zeros(convert(Int32, (landscape_size / dx) + 1), 
    convert(Int32, (landscape_size / dx) + 1));
  seed_disp_mat_builder_1D!(seed_disp_mat_1D, res = dx);
  #test to check the edge rows and cols sum to 0.5
  seed_disp_colsum = sum(seed_disp_mat_1D, 1);
  seed_disp_rowsum = sum(seed_disp_mat_1D, 2);
  
  pollen_disp_mat_1D = zeros(convert(Int32, (landscape_size / dx) + 1), 
    convert(Int32, (landscape_size / dx) + 1));
  pollen_disp_mat_builder_1D!(pollen_disp_mat_1D, res = dx);
  #test to check the edge rows and cols sum to 0.5
  pollen_disp_colsum = sum(pollen_disp_mat_1D, 1);
  pollen_disp_rowsum = sum(pollen_disp_mat_1D, 2);

#test the number of seeds dispersed to each location is the expected amount 

#test pollen dispersal and intergration works

#test number of survivours is as expected for a given worked example 

#test matrix mult gene mixing produces same result as the for loop version, 

Test.with_handler(cust_hand) do
  #test the disperal matrix is as expected 
  @test all(abs(seed_disp_colsum[[1 end]] - [0.5 0.5]) .<= 0.0000001)
  @test seed_disp_colsum[convert(Int32, floor(length(seed_disp_colsum) / 2))] >= 0.999
  @test all(abs(seed_disp_rowsum[[1 end]] - [0.5 0.5]) .<= 0.0000001)
  @test seed_disp_rowsum[convert(Int32, floor(length(seed_disp_colsum) / 2))] >= 0.999
  @test seed_disp_mat_1D[1, :] == transpose(seed_disp_mat_1D[:, 1])
  @test all(abs(pollen_disp_colsum[[1 end]] - [0.5 0.5]) .<= 0.0001)
  @test pollen_disp_colsum[convert(Int32, floor(length(pollen_disp_colsum) / 2))] >= 0.999
  @test all(abs(pollen_disp_rowsum[[1 end]] - [0.5 0.5]) .<= 0.0001)
  @test pollen_disp_rowsum[convert(Int32, floor(length(pollen_disp_colsum) / 2))] >= 0.999
  @test pollen_disp_mat_1D[1, :] == transpose(pollen_disp_mat_1D[:, 1])
  
  #same set of tests for the pollen kernel

end








#BENCH MARKING

#test the speed of matrix mult vs loop version of gene mixing

@time 