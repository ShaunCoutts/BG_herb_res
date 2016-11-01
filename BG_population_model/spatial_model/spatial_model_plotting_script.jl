# plotting functions for the spatial model and simulation experiments
# Will use the hiogh level interface Plots (which leans on multipule plotting backends
# this requires a plotting ecosystme to be imported, so a few comments to get the 
# packages I might need
# Pkg.add("Plots")
# Pkg.add("PyPlot")
# Pkg.add("GR")
# Pkg.add("Plotly") # get a build error on this one, don't load for now
# Pkg.add("StatPlots")
# Pkg.add("PlotRecipes")
# Pkg.add("Colors")

# load Plots, this should call the other backends as needed so 
# no need to load them
using Plots
using Colors
# choose a backend for Plots to wrap, I use PyPlot here for high
# levels of features 
pyplot()

# A function that takes the tuple of results output by 
# run_scen_trans, should make a 2 x 3 grid of plots
# one coloumn for each reciving landscape
# G_cols needs to be a row vector of 3 colours e.g. [:red :green :blue] corrosponding to [rr Rr RR]
function plot_scenario(output_tup::Tuple{Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}},
  G_cols, sink_col, burnin::Int64, output_loc::AbstractString, figure_name::AbstractString,
  adjust_scale::Float64)
  
  # go through the model output and turn any NaN cuased by the population not existing to a 0.0
  for i in 1:length(output_tup)
  
    output_tup[i][9, isnan(output_tup[i][7, :])] = NaN
    output_tup[i][5, output_tup[i][1, :] .== 0.0] = 0.0
    output_tup[i][6, output_tup[i][2, :] .== 0.0] = 0.0
    output_tup[i][7, output_tup[i][3, :] .== 0.0] = 0.0
    output_tup[i][8, output_tup[i][4, :] .== 0.0] = 0.0
  
  end
  
  # set up the fonts fo all the labels first 
  ax_font = Plots.Font("FreeSans", 12, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  title_font = Plots.Font("FreeSans", 14, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  leg_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  tic_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  
  # set up grib of 5 plots (to start, need to expand to 20 plots eventually 
  layout_arr = @layout [grid(1, 2); grid(1, 3)]
  plt = plot([0 0 0 0 0; burnin burnin burnin burnin burnin], ones(2, 5) * 1.1, #set up the color mat to show the burnin period
    fillrange = 0, fillalpha = 0.25, fillcolor = :lightgoldenrod, color = :white, label = "", 
    layout = layout_arr, xlim = (0, 70), ylim = (0, 1), grid = false, legend = :bottomright, #set up the layout and axis fonts
    title = ["a) Target site resistance" "b) Quantitative resistance" "c) Empty landscape" "d) Naive landscape" "e) Exposed landscape"], 
    titleloc = :left, xlabel = "time",
    ylabel = ["% R" "survival under g" "% occupied" "" ""],
    titlefont = title_font, guidefont = ax_font, tickfont = tic_font, 
    size = (600 * adjust_scale, 400 * adjust_scale))
    
  #start adding the other plots
  # % R plots
  pro_R = zeros(size(output_tup[1])[2], 3)
  pro_R[:, 1] = get_pro_R(output_tup[1][1, :], output_tup[1][2, :], output_tup[1][3, :]) #empty_%R
  pro_R[:, 2] = get_pro_R(output_tup[2][1, :], output_tup[2][2, :], output_tup[2][3, :]) #naive_%R
  pro_R[:, 3] = get_pro_R(output_tup[3][1, :], output_tup[3][2, :], output_tup[3][3, :]) #expos_%R
  
  ann_lab = Plots.PlotText("pre-TSR introduction", Plots.Font("FreeSans", 10, :hleft, :vcenter, 0.0, RGB{U8}(0.58, 0.644, 0.125)))
  
  plot!(plt, pro_R, subplot = 1, color = sink_col, label = ["empty" "naive" "exposed"], legend = :bottomright,
    annotations = (2, 0.9, ann_lab), legendfont = leg_font, linewidth = 3)
 
  # survival under g plots
  plot!(plt, [output_tup[1][9, :] output_tup[2][9, :] output_tup[3][9, :] ], subplot = 2, color = sink_col, legend = false, linewidth = 3)
    
  # %occ empty landscape
  plot!(plt, [output_tup[1][12, :] output_tup[1][11, :] output_tup[1][10, :] ], subplot = 3, color = G_col, legend = :topleft,
    label = ["rr" "Rr" "RR"], legendfont = leg_font, linewidth = 3)
    
  # %occ niave landscape 
  plot!(plt, [output_tup[2][12, :] output_tup[2][11, :] output_tup[2][10, :] ], subplot = 4, color = G_col, legend = false, linewidth = 3)
  
  # %occ expos landscape
  plot!(plt, [output_tup[3][12, :] output_tup[3][11, :] output_tup[3][10, :] ], subplot = 5, color = G_col, legend = false, linewidth = 3)
  
  # save the figure 
  #get present directory top return the working folder later
  int_file_loc = pwd()
  #set output locaiton
  cd(output_loc)
  savefig(plt, figure_name)
  cd(int_file_loc)
  
end

# function to produce a sinlge plot of %R over time and survival under mean g
function get_pro_R(RR_ts::Array{Float64, 1}, Rr_ts::Array{Float64, 1}, 
  rr_ts::Array{Float64, 1})
  
  return (2 * RR_ts + Rr_ts) ./ (2 * (RR_ts + Rr_ts + rr_ts))
  
end

function get_pro_R(RR::Float64, Rr::Float64, rr::Float64)
  
  return (2 * RR + Rr) / (2 * (RR + Rr + rr))
  
end

# function to produce a sinlge plot of %R over time and survival under mean g
function get_mean_spread(ts::Array{Float64, 1}, thresh::Float64)
  
  unsat_ts = ts[ts .< thresh]
  
  if(length(unsat_ts) > 2)
  
   return occ_dif = mean(unsat_ts[2:end] - unsat_ts[1:(end - 1)]) 
  
  else
  
    return(0)
  
  end
  
end

# produces a tuple that summaries the output_tup with a few snapshot numbers of the population in the final timestep
function trans_ex_snapshot(output_tup::Tuple{Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}},
  landscape_size::Float64, threshold::Float64)
  # set up a matrix to hold the results
  # source scenarios are in rows (empty, naive, exposed) and metrics are in coloumns
  # (%R, mean_spread_RRorRr, mean_spread_RR, mean_spread_Rr, mean_spread_rr) 
  snapshot = zeros(3, 5)
  
  # % R 
  snapshot[1, 1] = get_pro_R(output_tup[1][1, end], output_tup[1][2, end], output_tup[1][3, end]) #empty_%R
  snapshot[2, 1] = get_pro_R(output_tup[2][1, end], output_tup[2][2, end], output_tup[2][3, end]) #naive_%R
  snapshot[3, 1] = get_pro_R(output_tup[3][1, end], output_tup[3][2, end], output_tup[3][3, end]) #expos_%R
 
  # mean spread rates
  # reconstruct the number of locations occupied at each time step, rahter than proportion
  empty_occ = output_tup[1][10:12, burnin:end] * landscape_size 
  naive_occ = output_tup[2][10:12, burnin:end] * landscape_size 
  expos_occ = output_tup[3][10:12, burnin:end] * landscape_size 
  
  thres = threshold * landscape_size
  
  # do the RrorRR first, take max of Rr and RR at each time
  snapshot[1, 2] = get_mean_spread(max(empty_occ[1, :], empty_occ[2, :]), thres)
  snapshot[2, 2] = get_mean_spread(max(naive_occ[1, :], naive_occ[2, :]), thres)
  snapshot[3, 2] = get_mean_spread(max(expos_occ[1, :], expos_occ[2, :]), thres)
 
  # do the mean spread for each G
  for i in 1:3
    snapshot[1, i + 2] = get_mean_spread(empty_occ[i, :], thres)
    snapshot[2, i + 2] = get_mean_spread(naive_occ[i, :], thres)
    snapshot[3, i + 2] = get_mean_spread(expos_occ[i, :], thres)
  end
  
  # flatten out the snap shot matrix so it can fit in a row with the parameter values for easy plotting
  # order is [empty; naive; expos]
  return [snapshot[1, :]; snapshot[2, :]; snapshot[3, :]]      
  
end



