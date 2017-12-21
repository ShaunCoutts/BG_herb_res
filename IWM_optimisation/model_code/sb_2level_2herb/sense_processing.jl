# Proccess the sensitvity analysis results to a data frame
using HDF5, JLD

code_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/model_code/sb_2level_2herb"
data_loc = "/home/shauncoutts/Dropbox/projects/MHR_blackgrass/IWM_optimisation/outputs/GA_2level_sb_out/data_out"

cd(code_loc)
include("pop_process_TSR.jl"); 
include("act_seq_handelers.jl")
include("managment_functions.jl"); 
include("sense_helper_functions.jl")
include("plotting_functions_TSR.jl")

cd(data_loc)
sol_sweep = load("sens_runs_obj.jld")["sens_obj"];






