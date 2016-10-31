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





