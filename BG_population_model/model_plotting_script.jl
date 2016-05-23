# plotting the spatial model to see what is going on and start to understand what is happening 
using Colors
#using Gadfly
using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

using Cairo
using HDF5 JDL #for reading data objects from the disk

# the model output will be complex and large, having to be accesed from disk, so I define a 
# set of handler functions to extract the data 

# get mean g for each time (x) and locaiton (y), for each genotype  
function unpack_g(out_mat::Array{Float64,3})
  max_t = size(out_mat)[3]
  max_x = size(out_mat)[2]
  
  RR_g = reshape(out_mat[1, :, :], (max_x, max_t))
  Rr_g = reshape(out_mat[4, :, :], (max_x, max_t))
  rr_g = reshape(out_mat[7, :, :], (max_x, max_t))
  
  return (RR_g, Rr_g, rr_g)
  
end

function unpack_pop(out_mat::Array{Float64, 3})
  max_t = size(out_mat)[3]
  max_x = size(out_mat)[2]
  
  RR_pop = reshape(out_mat[3, :, :], (max_x, max_t))
  Rr_pop = reshape(out_mat[6, :, :], (max_x, max_t))
  rr_pop = reshape(out_mat[9, :, :], (max_x, max_t))
  
  return (RR_pop, Rr_pop, rr_pop)
 
end

# function to take the data object and produce a plot of three matricies and a legend
function spatial_plot(data_block::Tuple{Array{Float64, 1}, Array{Float64, 3}, Array{Float64, 3}}; 
  plot_herb::Bool = true)
  
  if plot_herb
    plot_mat = data_block[3]
  else
    plot_mat = data_block[2]
  end
  
  max_loc = size(plot_mat)[2]
  max_time = size(plot_mat)[3]
  loc_vals = 0:max_loc #number of coloumns in a slice
  time_vals = 0:max_time #number of time slices
  min_g = minimum(plot_mat[[1, 4, 7], :, :])
  max_g = maximum(plot_mat[[1, 4, 7], :, :])
  max_pop = maximum(plot_mat[[3, 6, 9], :, :])
  # set some figure margin constants
  mar_left = 0.1
  mar_width = 0.8
  mar_height = 0.2
  mar_bots = [0.78, 0.56, 0.34, 0.06]
 
  #get some nicer tuples to work with 
  g_tup = unpack_g(plot_mat)
  pop_tup = unpack_pop(plot_mat)
  
  #build figure
  fig_ob = figure("MR_TSR_space_plot", figsize = (10, 10))
  
  # set up the RR subplot
  ax_RR = fig_ob[:add_subplot](4, 1, 1)
  new_pos = [0.06, 0.7, 0.70, 0.9]
  new_pos = [mar_left, mar_bots[1], mar_width, mar_height]
  ax_RR[:set_position](new_pos)
  setp(ax_RR[:get_xticklabels](),visible=false) # Disable x tick labels
  setp(ax_RR[:get_yticklabels](),visible=false) # Disable x tick labels
  ylabel("RR\nlocation")
  ylim(0.0, max_loc)
  xlim(0.0, max_time)
  yticks(0 : 50 : max_loc)
  
  # set up the Rr subplot
  ax_Rr = fig_ob[:add_subplot](4, 1, 2)
  new_pos = [mar_left, mar_bots[2], mar_width, mar_height]
  ax_Rr[:set_position](new_pos)
  setp(ax_Rr[:get_xticklabels](),visible=false) # Disable x tick labels
  setp(ax_Rr[:get_yticklabels](),visible=false) # Disable x tick labels
  ylabel("Rr\nlocation")
  ylim(0.0, max_loc)
  xlim(0.0, max_time)
  yticks(0 : 50 : max_loc)
 
  #set up the rr subplot
  ax_rr = fig_ob[:add_subplot](4, 1, 3)
  new_pos = [mar_left, mar_bots[3], mar_width, mar_height]
  ax_rr[:set_position](new_pos)
  setp(ax_rr[:get_yticklabels](),visible=false) # Disable x tick labels
  ylabel("rr\nlocation")
  ylim(0.0, max_loc)
  xlim(0.0, max_time)
  yticks(0 : 50 : max_loc)
  xlabel("time step")
  
  #set up the axis, ideally a square 
  ax_leg = fig_ob[:add_subplot](4, 1, 4)
  new_pos = [mar_left, mar_bots[4], 0.2, 0.2]
  ax_leg[:set_position](new_pos)
  ylabel("res. score (g)")
  xlabel("density (max = 1)")
  yticks(min_g : ((max_g - min_g) / 4) : max_g)
  xticks(0 : 0.25 : 1)
  title("color map legend")
  #create a mat of g values and pop values to plot colours for
  dummy_g_mat = hcat([collect(linspace(min_g, max_g, 100)) for i in 1:100]...)
  dummy_pop_mat = transpose(hcat([collect(linspace(0, max_pop, 100)) for i in 1:100]...))
  plot_col_mat!(ax_leg, dummy_g_mat, dummy_pop_mat; min_val_1 = min_g, max_val_1 = max_g, 
    min_val_2 = 0.0, max_val_2 = max_pop, min_x = 0.0, max_x = 1.0, min_y = min_g, max_y = max_g)
  
end

# home spun function that takes two values and returns a HSL colour where val_1
# controls hue and val_2 controls chroma
function two_channle_col(val_1::Float64, val_2::Float64, min_val_1::Float64 = 0.0,
  max_val_1::Float64 = 1.0, min_val_2::Float64 = 0.0, max_val_2::Float64 = 1.0)
  
  rescale_hue = (270.0 * (val_1 - min_val_1)) / (max_val_1 - min_val_1)
  
  min_light, max_light = 0.15, 1 #0 = all black, 1 = all white
  #rescale and flip so smaller values get whiter
  rescale_light = (((max_light - min_light) * (max_val_2 - val_2)) / (max_val_2 - min_val_2)) + min_light
  
  col = convert(RGB, HSL(rescale_hue, 1.0, rescale_light))
  
  return col
  
end

# plot the color matrix and add it to an existing figure object
# function that takes a two matrices, of numbers of var1 (controls colors red to bule)
# and a second matrix of var2 that controls the croma (amount of color, from white to 
# fully colored). Return a matrix of colors
function plot_col_mat!(ax_ob::PyCall.PyObject, val_mat_1::Array{Float64, 2}, val_mat_2::Array{Float64, 2}; 
  min_val_1::Float64 = 0.0, max_val_1::Float64 = 1.0, min_val_2::Float64 = 0.0, max_val_2::Float64 = 1.0,
  min_x::Float64 = 0.0, max_x::Float64 = 1.0, min_y::Float64 = 0.0, max_y::Float64 = 1.0)

  seg_height = 1 / size(val_mat_1)[1]
  seg_width = 1 / size(val_mat_1)[2]
  
  x_coord = (((max_x - min_x) * (collect(1:size(val_mat_1)[2]) - 1)) / (size(val_mat_1)[2] - 1)) + min_x
  y_coord = (((max_x - min_x) * (collect(1:size(val_mat_1)[1]) - 1)) / (size(val_mat_1)[1] - 1)) + min_x
  
  for x in 1:size(val_mat_1)[2]
    for y in 1:size(val_mat_1)[1]
#       col = convert(RGB, HSL(0.0, 0.5, 0.5)) 
#       col_tup = (red(col), green(col), blue(col))
#      
      #rescale the matrix space to the co-ornate syastme of the plot
      x_pos = ()
      col = two_channle_col(val_mat_1[y, x], val_mat_2[y, x], min_val_1, max_val_1,
	min_val_2, max_val_2)
      col_tup = (red(col), green(col), blue(col))
      rect = patch.Rectangle((x_coord[x], y_coord[y]), seg_width, seg_height, linewidth = 0.0, facecolor = col_tup)
      ax_ob[:add_artist](rect)
    end
  end
  
  return nothing
 
end


