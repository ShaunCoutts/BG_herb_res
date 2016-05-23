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
 
  fig_ob = figure("MR_TSR_space_plot", figsize = (10, 10))
#   subplots_adjust(hspace = 0.2)
  
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
  
#   two_channle_mat()
#   plot_col_mat!(ax)#, col_mat)
#   
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
  
end

# function that takes a two matrices, of numbers of var1 (controls colors red to bule)
# and a second matrix of var2 that controls the croma (amount of color, from white to 
# fully colored). Return a matrix of colors
function two_channle_mat(var1::Array{Float64, 2}, var2::Array{Float64, 2})

   
end

# plot the color matrix and add it to an existing figure object
function plot_col_mat!(ax_ob::PyCall.PyObject)#, col_mat::Array{Colours, 2})

  seg_height = 0.1
  seg_width = 0.1
  
  xy_co = ([0.1, 0.1], [0.1, 0.3])

 
  for corn in xy_co
      rect = patch.Rectangle((corn[1], corn[2]), seg_width, seg_height)#, fill = col_mat)
      ax_ob[:add_artist](rect)
  end
  
  return nothing
 
end


