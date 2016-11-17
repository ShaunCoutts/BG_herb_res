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
using StatPlots
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
  
  return nothing
  
end


# visulisation of univiriant parameter sweep to see how parameter vlaues affect
# %R and survival under g at final time step, first attempt showing 6 lines on a single 
# plot, this may be a bit much
function pop_res_4_scen(res_df::DataFrame, x_var::Symbol, y_var::Array{Symbol, 1},  
  sink_col, output_loc::AbstractString, figure_name::AbstractString, adjust_scale::Float64, 
  x_lab::String, y_lab::Array{String, 1}, low_TSR_inj::Float64, low_g_inj::Float64, high_TSR_inj::Float64, 
  high_g_inj::Float64)
 
  # get a little basic info about the plotting axes
  x_min = floor(minimum(res_df[x_var]), 1)
  x_max = ceil(maximum(res_df[x_var]), 1)
  y_min = floor(minimum(vcat(colwise(minimum, res_df[y_var]) ...)), 1)
  y_max = ceil(maximum(vcat(colwise(maximum, res_df[y_var]) ...)), 1)
  
  #rearrange the sink colours to reflect the order they are plotted from the dataframe
  sink_col_ra = sink_col[[1 3 2]]
  #make a grey scale to choose some shades from 
  grey_pal = colormap("Grays", 100)
  
  # set up the fonts fo all the labels first 
  ax_font = Plots.Font("FreeSans", 12, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  title_font = Plots.Font("FreeSans", 14, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  leg_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  tic_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  
  # make a grid of 8 empty plots with all the formatting on them 
  layout_arr = @layout [grid(2, 1) grid(2, 1); grid(2, 1) grid(2, 1)]
  plt = plot(xlim = (x_min, x_max), ylim = (y_min, y_max), layout = layout_arr, grid = false, legend = :bottomleft,
    title = ["a) quant. res. low, TSR low" "" "b) quant. res. low. TSR high" "" "c) quant. res. high, TSR low" "" "d) quant. res. high, TSR high" ""],
    titleloc = :left, titlefont = title_font, guidefont = ax_font, tickfont = tic_font, legendfont = leg_font, 
    size = (600 * adjust_scale, 400 * adjust_scale), border = false, bg_inside = grey_pal[10])
  
  # plot the low g low TSR scen
  plot!(plt, res_df[(res_df[:inj_TSR] .== low_TSR_inj) & (res_df[:inj_g] .== low_g_inj), :], x_var, y_var[1], group = :scen, 
    subplot = 1, color = sink_col_ra, linewidth = 3, xlabel = "", ylabel = y_lab[1], 
    label = ["empty-%R" "expos-%R" "naive-%R"])

  if length(y_var) == 2
    plot!(plt, res_df[(res_df[:inj_TSR] .== low_TSR_inj) & (res_df[:inj_g] .== low_g_inj), :], x_var, y_var[2], group = :scen, 
      subplot = 2, color = sink_col_ra, linewidth = 3, linestyle = :dash, 
      label = ["empty-sur rr" "expos-sur rr" "naive-sur rr"], ylabel = y_lab[2], xlabel = "")
  end
  
 
  # plot the low g, high TSR scen 
  plot!(plt, res_df[(res_df[:inj_TSR] .== high_TSR_inj) & (res_df[:inj_g] .== low_g_inj), :], x_var, y_var[1], group = :scen, 
    subplot = 3, color = sink_col_ra, linewidth = 3, xlabel = "", ylabel = "", legend = false)
 
  if length(y_var) == 2
    plot!(plt, res_df[(res_df[:inj_TSR] .== high_TSR_inj) & (res_df[:inj_g] .== low_g_inj), :], x_var, y_var[2], group = :scen, 
      subplot = 4, color = sink_col_ra, linewidth = 3, xlabel = "", ylabel = "",  linestyle = :dash, legend = false)
  end
    
   # plot the high g, low TSR scen 
  plot!(plt, res_df[(res_df[:inj_TSR] .== low_TSR_inj) & (res_df[:inj_g] .== high_g_inj), :], x_var, y_var[1], group = :scen, 
    subplot = 5, color = sink_col_ra, linewidth = 3, xlabel = "", ylabel = y_lab[1], legend = false)
    
  if length(y_var) == 2
    plot!(plt, res_df[(res_df[:inj_TSR] .== low_TSR_inj) & (res_df[:inj_g] .== high_g_inj), :], x_var, y_var[2], group = :scen, 
      subplot = 6, color = sink_col_ra, linewidth = 3, linestyle = :dash, legend = false, xlabel = x_lab, ylabel = y_lab[2])
  end
   
  # plot the high g, high TSR scen 
 
  plot!(plt, res_df[(res_df[:inj_TSR] .== high_TSR_inj) & (res_df[:inj_g] .== high_g_inj), :], x_var, y_var[1], group = :scen, 
    subplot = 7, color = sink_col_ra, linewidth = 3, xlabel = "", ylabel = "", legend = false)
  if length(y_var) == 2
    plot!(plt, res_df[(res_df[:inj_TSR] .== high_TSR_inj) & (res_df[:inj_g] .== high_g_inj), :], x_var, y_var[2], group = :scen, 
      subplot = 8, color = sink_col_ra, linewidth = 3, linestyle = :dash, legend = false, xlabel = x_lab, ylabel = "")
  end
    
  # save the figure 
  #get present directory top return the working folder later
  int_file_loc = pwd()
  #set output locaiton
  cd(output_loc)
  savefig(plt, figure_name)
  cd(int_file_loc)
  
  return nothing
    
end

# make a plot for the change in TSR injection number 
function TSR_inj_var_plot(res_df::DataFrame, x_var::Symbol, y_var::Array{Symbol, 1},  
  sink_col, output_loc::AbstractString, figure_name::AbstractString, adjust_scale::Float64, 
  x_lab::String, y_lab::Array{String, 1}, low_g_inj::Float64, high_g_inj::Float64)
  
  # get a little basic info about the plotting axes
  x_min = floor(minimum(res_df[x_var]), 1)
  x_max = ceil(maximum(res_df[x_var]), 1)
  y_min = floor(minimum(vcat(colwise(minimum, res_df[y_var]) ...)), 1)
  y_max = ceil(maximum(vcat(colwise(maximum, res_df[y_var]) ...)), 1)
  
  #rearrange the sink colours to reflect the order they are plotted from the dataframe
  sink_col_ra = sink_col[[1 3 2]]
  #make a grey scale to choose some shades from 
  grey_pal = colormap("Grays", 100)
  
  # set up the fonts fo all the labels first 
  ax_font = Plots.Font("FreeSans", 12, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  title_font = Plots.Font("FreeSans", 14, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  leg_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  tic_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  
  # make a grid of 4 empty plots with all the formatting on them 
  layout_arr = @layout grid(2, 2)
  plt = plot(xlim = (x_min, x_max), ylim = (y_min, y_max), layout = layout_arr, grid = false,
    title = ["a) quant. res. low" "b) quant. res. high" "" ""],
    titleloc = :left, titlefont = title_font, guidefont = ax_font, tickfont = tic_font, legendfont = leg_font, 
    size = (600 * adjust_scale, 600 * adjust_scale), border = false, bg_inside = grey_pal[10])
 
  # plot the low g
  plot!(plt, res_df[res_df[:inj_g] .== low_g_inj, :], x_var, y_var[1], group = :scen, 
    subplot = 1, color = sink_col_ra, linewidth = 3, xlabel = "", ylabel = y_lab[1], legend = :bottomleft)
 
  plot!(plt, res_df[res_df[:inj_g] .== low_g_inj, :], x_var, y_var[2], group = :scen, 
    subplot = 3, color = sink_col_ra, linewidth = 3, xlabel = x_lab, ylabel = y_lab[2], linestyle = :dash, 
    legend = :topright)
   
  # plot the low g
  plot!(plt, res_df[res_df[:inj_g] .== high_g_inj, :], x_var, y_var[1], group = :scen, 
    subplot = 2, color = sink_col_ra, linewidth = 3, xlabel = "", ylabel = "", legend = false)
 
  plot!(plt, res_df[res_df[:inj_g] .== high_g_inj, :], x_var, y_var[2], group = :scen, 
    subplot = 4, color = sink_col_ra, linewidth = 3, xlabel = x_lab, ylabel = "",  linestyle = :dash, 
    legend = false)
  
  # save the figure 
  #get present directory top return the working folder later
  int_file_loc = pwd()
  #set output locaiton
  cd(output_loc)
  savefig(plt, figure_name)
  cd(int_file_loc)
  
  return nothing

end
  
# produce 3 2D colormat plots to show how 2 varaibles (:x, :y) interact to determine the 3rd (:z).
function colormats_2D(res_df::DataFrame, x_var::Symbol, y_var::Symbol, z_var::Symbol,  
  output_loc::AbstractString, figure_name::AbstractString, adjust_scale::Float64, 
  x_lab::String, y_lab::String, z_lab::String, inj_scen_name::String)
  
  #make a grey scale to choose some shades from 
  grey_pal = colormap("Grays", 100)
  
  # set up the fonts fo all the labels first 
  ax_font = Plots.Font("FreeSans", 12, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  title_font = Plots.Font("FreeSans", 14, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  leg_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  tic_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0))
  
  # make a grid of 4 empty plots with all the formatting on them 
  layout_arr = @layout [b{0.0001h}; a{0.25w} a{0.5w} a{0.25w}; grid(1, 2)]
  plt = plot(layout = layout_arr, grid = false, background_color_outside = grey_pal[10], border = false,
    background_color = grey_pal[10], title = [inj_scen_name "" "a) empty" "" "b) naive" "c) exposed"],
    titleloc = :left, titlefont = title_font, guidefont = ax_font, tickfont = tic_font, legendfont = leg_font, 
    size = (700 * adjust_scale, 600 * adjust_scale), xlabel = ["" "" x_lab "" x_lab x_lab], 
    ylabel = ["" "" y_lab "" y_lab ""])
 
  # set up the plotting axes
  x_ax = unique(res_df[x_var])
  y_ax = unique(res_df[y_var])
  
  # make a contour plot for each reciving plot_scenario
  df_empty = convert(Array{Any, 2}, res_df[res_df[:scen] .== "empty", [x_var, y_var, z_var]])
  df_naive = convert(Array{Any, 2}, res_df[res_df[:scen] .== "naive", [x_var, y_var, z_var]])
  df_expos = convert(Array{Any, 2}, res_df[res_df[:scen] .== "expos", [x_var, y_var, z_var]])
  
  z_mat_empty = zeros(length(y_ax), length(x_ax))
  z_mat_naive = zeros(length(y_ax), length(x_ax))
  z_mat_expos = zeros(length(y_ax), length(x_ax))
  
  for x in 1:length(x_ax)
    for y in 1:length(y_ax)
      
      z_mat_empty[y, x] = z_finder(x_ax[x], y_ax[y], df_empty)[1] 
      z_mat_naive[y, x] = z_finder(x_ax[x], y_ax[y], df_naive)[1] 
      z_mat_expos[y, x] = z_finder(x_ax[x], y_ax[y], df_expos)[1] 
      
    end
  end
 
  contour!(x_ax, y_ax, z_mat_empty, fill = true, subplot = 3) 
  contour!(x_ax, y_ax, z_mat_naive, fill = true, subplot = 5) 
  contour!(x_ax, y_ax, z_mat_expos, fill = true, subplot = 6) 
  
  # save the figure 
  #get present directory top return the working folder later
  int_file_loc = pwd()
  #set output locaiton
  cd(output_loc)
  savefig(plt, figure_name)
  cd(int_file_loc)
  
  return nothing

end

# function to take two values and find a thrid value to return, 
# order of search array should be [x y z]
function z_finder(targ_x, targ_y, search_array::Array{Any, 2})
  
  return search_array[(search_array[:, 1] .== targ_x) & (search_array[:, 2] .== targ_y), 3]
  
end  

function logit_sur_2_prob(logit_vect, base_sur::Float64)

  return 1 ./ (1 + exp(-(base_sur - logit_vect))) 

end

##############################################################################################
## PLOTTING FOR THE NATURAL SPREAD EXPERIMENTS ###############################################

# take 2 matricies of floats, produce a matrix of HSL colors using 2 channles, 
# the first channle controls lightness, 2nd channle does hue 
function colmat_2channel(mat1::Array{Float64, 2}, mat2::Array{Float64, 2}, 
  min_light::Float64, max_light::Float64, hue_start::Float64, hue_end::Float64, 
  min_c1::Float64, max_c1::Float64, min_c2::Float64, max_c2::Float64)
  
  # map all vlaues in mat1 to scale between min_light and max_light
  mat_dim = size(mat1) # use same size for both mats so it throws an error if they are different shapes
  col_mat = Array{RGB{Float64}, 2}(mat_dim)
  for y in 1:mat_dim[1]
    for x in 1:mat_dim[2]
      
      light = rescale(mat1[y, x], min_c1, max_c1, min_light, max_light)
      hue = rescale(mat2[y, x], min_c2, max_c2, hue_start, hue_end)
      col_mat[y, x] = convert(RGB, HSL(hue, 1.0, light))
   
    end
  end
  
  return col_mat
  
end

function rescale(x::Float64, old_min::Float64, old_max::Float64,
  new_min::Float64, new_max::Float64)
  
  return (((new_max - new_min) * (x - old_min)) / (old_max - old_min)) + new_min  
  
end

function legend_colmat(min_x::Float64, max_x::Float64, min_y::Float64, max_y::Float64,
  min_light::Float64, max_light::Float64, hue_start::Float64, hue_end::Float64)
  
  resolution = 100
  x_ax = transpose(repeat(linspace(min_x, max_x, resolution), outer = (1, resolution)))
  y_ax = repeat(linspace(min_y, max_y, resolution), outer = (1, resolution))
  y_ax = y_ax[end:-1:1, :]
  
  col_leg = colmat_2channel(x_ax, y_ax, min_light, max_light, hue_start, hue_end, 
    min_x, max_x, min_y, max_y)
  
 return col_leg 
  
end

# grid of two channel heat maps to show how populations change over time
# pass in 6 colormatircies to plot 
function dualchan_heatmap_grid(ee_empty::Array{RGB{Float64}, 2}, ee_full::Array{RGB{Float64}, 2},
  ne_empty::Array{RGB{Float64}, 2}, ne_full::Array{RGB{Float64}, 2}, en_empty::Array{RGB{Float64}, 2},
  en_full::Array{RGB{Float64}, 2}, adjust_scale::Float64, light_max::Float64, z_min::Float64, z_max::Float64, 
  z_lab::String, hue_min::Float64, hue_max::Float64, output_loc::String, output_name::String)

  #make a grey scale to choose some shades from 
  grey_pal = colormap("Grays", 100);

  # set up the fonts fo all the labels first 
  ax_font = Plots.Font("FreeSans", 12, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0));
  title_font = Plots.Font("FreeSans", 14, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0));
  leg_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0));
  tic_font = Plots.Font("FreeSans", 10, :hcenter, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0));
  
  # bit of annotation to label regions
  ann_lab_source = Plots.PlotText("source of TSR", Plots.Font("FreeSans", 10, :hleft, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0)))
  ann_lab_recive = Plots.PlotText("reciving area", Plots.Font("FreeSans", 10, :hleft, :vcenter, 0.0, RGB{U8}(0.0, 0.0, 0.0)))

  # now make the plots 
  layout_arr = @layout [grid(2, 3) a{0.1w}]; 
    
  plt = plot(layout = layout_arr, grid = false, background_color_outside = grey_pal[10], 
    border = false, background_color = grey_pal[10],   
    title = ["a) source exposed\n    recive exposed" "b) source naive\n    recive exposed" "c) source exposed\n    recive naive"  "d)" "e)" "f)" ""],
    titleloc = :left, titlefont = title_font, guidefont = ax_font, tickfont = tic_font, 
    legendfont = leg_font, size = (600 * adjust_scale, 400 * adjust_scale), 
    xlabel = ["" "" "" "time (years)" "time (years)" "tiome (years)" ""], 
    ylabel = ["space (m)" "" "" "space (m)" "" "" z_lab]);

  heatmap!(plt, ee_empty, subplot = 1, annotations = [(10, 10, ann_lab_source), (10, 30, ann_lab_recive)], 
    xticks = []);
  heatmap!(plt, ee_full, subplot = 4);

  heatmap!(plt, ne_empty, subplot = 2, xticks = [], yticks = []);
  heatmap!(plt, ne_full, subplot = 5, yticks = []);

  heatmap!(plt, en_empty, subplot = 3, xticks = [], yticks = []);
  heatmap!(plt, en_full, subplot = 6, yticks = []);

  leg_mat = legend_colmat(light_max, light_max + 1, z_min, z_max, 0.5, 0.5, hue_min, hue_max)
  heatmap!(plt, leg_mat, subplot = 7, yticks = ([100, 75, 50, 25, 1], [z_min, (z_max * 0.25) + z_min, 
    (z_max * 0.5) + z_min, (z_max * 0.75) + z_min, z_max]), xticks = [], aspect_ratio = 3.0);
  #add line to show the source and sink
  for i in 1:6
    plot!(plt, [1, num_iter], [source_locs[end], source_locs[end]], subplot = i, label = "",
      linecolor = :black);
  end
  plot!(plt, [1, num_iter], [source_locs[end], source_locs[end]], subplot = 1, label = "",
      linecolor = :black);

  # save the figure 
  #get present directory top return the working folder later
  int_file_loc = pwd()
  #set output locaiton
  cd(output_loc)
  savefig(plt, output_name)
  cd(int_file_loc)

  return nothing

end
