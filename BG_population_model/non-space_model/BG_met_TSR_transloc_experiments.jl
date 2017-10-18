# simulation experiments for the TSR and NTSR joint evolution
# to test how TSR can invade into a population with NTSR
@everywhere using DistributedArrays

using DataFrames

# need to run each source scenario three times, once in a empty landscape, one in an exposed population
# and one in a herbicide exposed population (can pre calculate all these so don't need to re-caclulate each time 
@everywhere file_loc_func_p = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/non-space_model" 
output_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_population_model/model_output" 

@everywhere cd(file_loc_func_p)
@everywhere include("BG_met_TSR_pop_process.jl")
@everywhere include("BG_met_TSR_runners.jl")
# include("spatial_model_plotting_script.jl")
# @everywhere include("BG_met_TSR_exper_funs.jl")
  

###################################################################################################
# script to run the translocation experiments
# parameter values
upper_g = 20.0;
lower_g = -20.0;
dg = 0.5;
int_mean_g = 0.0;
int_sd_g = 1.4142;

int_num_RR = 0.0;
int_num_Rr = 0.0;
int_num_rr = 10.0; # number of intial seeds at each location for each genoptype, assume only TS susceptible
num_iter = 100;
int_run = 20;

offspring_sd = 1.0;
fec_max = 60.0;
fec0 = 4.0;
fec_cost = 0.45;
dd_fec = 0.001;
base_sur = 10.0; 
herb_effect = 16.0; 
g_prot = 1.5; 
pro_exposed = 0.8;
seed_sur = 0.45;
germ_prob = 0.52;
resist_G = ["RR", "Rr"];
herb_app = 2;

# set up the evaluation points for quantitative resistance
g_vals = collect(lower_g : dg : upper_g);   

# set a up a limited parameter to find a nice region of parameter space to work in
g_pro = [1.0, 1.5, 2.0];
Va = [0.5, 1.0, 1.5]; # addative variance, take sqrt() to standard deviation, which is the julia parameterisation

par_list = [];
rep_count = 1.0;
for gp in g_pro
  for osv in Va

    push!(par_list, [int_num_RR, int_num_Rr, int_num_rr, germ_prob, fec0, fec_cost, fec_max, 
      dd_fec, herb_effect, gp, seed_sur, pro_exposed, base_sur, sqrt(osv), rep_count]);
    
    rep_count += 1;
    
  end
end

# copy that list to all the  workers 
@eval @everywhere par_list = $par_list

# define the other needed inputs on all workers
@everywhere upper_g = 20.0;
@everywhere lower_g = -20.0;
@everywhere dg = 0.5;
@everywhere int_mean_g = 0.0;
@everywhere int_sd_g = 1.4142;
@everywhere num_iter = 100;
@everywhere resist_G = ["RR", "Rr"];
@everywhere herb_app = 2;
@everywhere g_vals = collect(lower_g : dg : upper_g);   


@time out = @DArray [run_wrapper(par_list[x], int_mean_g, int_sd_g, num_iter, g_vals, dg, 
  resist_G, herb_app, "no TSR") for x = 1:length(par_list)];


var_names = ["int_mean_g", "int_sd_g", "intRR", "intRr", "intrr", "germ_prob", "fec0", "fec_cost", 
  "fec_max", "dd_fec", "herb_effect", "g_pro", "seed_sur", "pro_exposed", "s0", "off_sd", "scen", "rep_ID", "measure"];
all_names = vcat(var_names, [string("t", i) for i = 1:num_iter]);

res_df = DataFrame(vcat(out...));
names!(res_df, convert(Array{Symbol}, all_names));

# write the parameter sweep to a .csv file so an R plotting script can be used
cd(output_loc);
writetable("rho_Va_par_sweep_nonspace.csv", res_df);
  
###############################################################################################################################################
## NOW RUN WITH TSR   
# define the other needed inputs on all workers
upper_g = 20.0;
lower_g = -20.0;
dg = 0.5;
num_est = 20;
num_iter = 100;
resist_G = ["RR", "Rr"];
g_vals = collect(lower_g : dg : upper_g);   

fec_max = 60.0;
fec0 = 4.0;
fec_cost = 0.45;
dd_fec = 0.0005;
base_sur = 10.0; 
herb_effect = 16.0; 
g_prot = collect(1.0 : 0.1 : 2.0); 
pro_exposed = 0.8;
seed_sur = 0.45;
germ_prob = 0.52;
resist_G = ["RR", "Rr"];

Va = [0.5, 1.0, 1.5]; # addative variance, take sqrt() to standard deviation, which is the julia parameterisation
herb_app1 = [1, 1, 2];
herb_app2 = 2;

int_num_RR = 0.0;
int_num_Rr = 0.0;
int_num_rr = [0.0, 1000.0, 1000.0];
int_g = [0.0, 0.0, get_g_at_sur(0.5, base_sur, herb_effect, g_prot[1])];
int_sd = [1.4142, 1.4142, 1.0];

inj_num_RR = [0.01, 0.1, 1.0, 10.0];
inj_num_Rr = 0.0;
inj_num_rr = 10.0 - inj_num_RR;
inj_sd = [1.4142, 1.0, 1.0];
s_sur_rr = [0.01, 0.5, 0.97]; # target survival for source population   
  
@time out = [run_wrapper_hot_seed_injection(int_num_RR, int_num_Rr, int_num_rr[tscen], 
  inj_num_RR[inj_R], inj_num_Rr, inj_num_rr[inj_R], int_g[tscen], int_sd[tscen], 
  get_g_at_sur(s_sur_rr[sg], base_sur, herb_effect, g_prot[rho]), inj_sd[sg], num_est, num_iter, 
  herb_app1[tscen], herb_app2, germ_prob, fec0, fec_cost, fec_max, dd_fec, herb_effect, g_prot[rho], 
  seed_sur, pro_exposed, base_sur, sqrt(Va[v]), g_vals, dg, resist_G) for 
  v = 1:length(Va), rho = 1:length(g_prot), tscen = 1:length(herb_app1), 
  inj_R = 1:length(inj_num_RR), sg = 1:length(s_sur_rr)];
  
var_names = ["intRR", "intRr", "intrr", "injRR", "injRr", "injrr", "int_g", "int_sd", "inj_g", "inj_sd", "herb1", "herb2",
  "germ_prob", "fec0", "fec_cost", "fec_max", "dd_fec", "herb_effect", "g_pro", "seed_sur", "pro_exposed", "s0", 
  "off_sd", "est_period", "measure"];
all_names = vcat(var_names, [string("t", i) for i = 1:(num_est + num_iter)]);

res_df = DataFrame(vcat(out...));
names!(res_df, convert(Array{Symbol}, all_names));

# write the parameter sweep to a .csv file so an R plotting script can be used
cd(output_loc);
writetable("rho_Va_trans_expr.csv", res_df);
  
## NO TSR RUN  
int_num_RR = 0.0;
int_num_Rr = 0.0;
int_num_rr = 0.0;
int_g = 0.0;
int_sd = 1.4142;

g_prot = [1.0, 1.2, 1.5]; 
Va = [0.5, 1.0, 1.5]; # addative variance, take sqrt() to standard deviation, which is the julia parameterisation
herb_app1 = 1;

inj_num_RR = 0.0;
inj_num_Rr = 0.0;
inj_num_rr = 10.0;
inj_sd = [1.4142, 1.0, 1.0];
s_sur_rr = [0.01, 0.5, 0.97]; # target survival for source population   
  
@time out = [run_wrapper_hot_seed_injection(int_num_RR, int_num_Rr, int_num_rr, 
  inj_num_RR, inj_num_Rr, inj_num_rr, int_g, int_sd, 
  get_g_at_sur(s_sur_rr[sg], base_sur, herb_effect, g_prot[rho]), inj_sd[sg], num_est, num_iter, 
  herb_app1, herb_app2, germ_prob, fec0, fec_cost, fec_max, dd_fec, herb_effect, g_prot[rho], 
  seed_sur, pro_exposed, base_sur, sqrt(Va[v]), g_vals, dg, resist_G) for 
  v = 1:length(Va), rho = 1:length(g_prot), sg = 1:length(s_sur_rr)];
  
var_names = ["intRR", "intRr", "intrr", "injRR", "injRr", "injrr", "int_g", "int_sd", "inj_g", "inj_sd", "herb1", "herb2",
  "germ_prob", "fec0", "fec_cost", "fec_max", "dd_fec", "herb_effect", "g_pro", "seed_sur", "pro_exposed", "s0", 
  "off_sd", "est_period", "measure"];
all_names = vcat(var_names, [string("t", i) for i = 1:(num_est + num_iter)]);

res_df = DataFrame(vcat(out...));
names!(res_df, convert(Array{Symbol}, all_names));

# write the parameter sweep to a .csv file so an R plotting script can be used
cd(output_loc);
writetable("rho_Va_noTSR.csv", res_df);
  

##########################################################################################################################################################################
# run to get the effect on TSR fitness advantage for the non-spatial model, for just 1 parameter set
upper_g = 20.0;
lower_g = -20.0;
dg = 0.5;
num_est = 20;
num_iter = 100;
resist_G = ["RR", "Rr"];
g_vals = collect(lower_g : dg : upper_g);   

fec_max = 60.0;
fec0 = 4.0;
fec_cost = 0.45;
dd_fec = 0.0005;
base_sur = 10.0; 
herb_effect = 16.0; 
g_prot = 1.5;
pro_exposed = 0.8;
seed_sur = 0.45;
germ_prob = 0.52;
resist_G = ["RR", "Rr"];

Va = [0.5, 1.0, 1.5]; # addative variance, take sqrt() to standard deviation, which is the julia parameterisation
herb_app1 = [1, 1, 2];
herb_app2 = 2;

int_num_RR = 0.0;
int_num_Rr = 0.0;
int_num_rr = [0.0, 1000.0, 1000.0];
int_g = [0.0, 0.0, get_g_at_sur(0.5, base_sur, herb_effect, g_prot)];
int_sd = [1.4142, 1.4142, 1.0];

inj_num_RR = 1.0;
inj_num_Rr = 0.0;
inj_num_rr = 10.0 - inj_num_RR;
inj_sd = 1.0;
s_sur_rr = 0.5; # target survival for source population   
  
@time out = [run_wrapper_hot_seed_injection(int_num_RR, int_num_Rr, int_num_rr[tscen], 
  inj_num_RR, inj_num_Rr, inj_num_rr, int_g[tscen], int_sd[tscen], 
  get_g_at_sur(s_sur_rr, base_sur, herb_effect, g_prot), inj_sd, num_est, num_iter, 
  herb_app1[tscen], herb_app2, germ_prob, fec0, fec_cost, fec_max, dd_fec, herb_effect, g_prot, 
  seed_sur, pro_exposed, base_sur, sqrt(v), g_vals, dg, resist_G) for 
  v = Va, tscen = 1:length(herb_app1)];
  
var_names = ["intRR", "intRr", "intrr", "injRR", "injRr", "injrr", "int_g", "int_sd", "inj_g", "inj_sd", "herb1", "herb2",
  "germ_prob", "fec0", "fec_cost", "fec_max", "dd_fec", "herb_effect", "g_pro", "seed_sur", "pro_exposed", "s0", 
  "off_sd", "est_period", "measure"];
all_names = vcat(var_names, [string("t", i) for i = 1:(num_est + num_iter)]);

res_df = DataFrame(vcat(out...));
names!(res_df, convert(Array{Symbol}, all_names));

# write the parameter sweep to a .csv file so an R plotting script can be used
cd(output_loc);
writetable("TSR_adv_trans_expr.csv", res_df);
  