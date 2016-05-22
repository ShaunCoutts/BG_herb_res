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

# function to take the data object and produce a plot
function spatial_plot(data_bloack::Tuple{Array{Float64, 1}, Array{Float64, 3}, Array{Float64, 3}}; 
  plot_herb::Bool = true)
  
  #plot stuff in here
  
end

# function that takes a two matrices, of numbers of var1 (controls colors red to bule)
# and a second matrix of var2 that controls the croma (amount of color, from white to 
# fully colored)

# or it just takes a matrix of colors and plots them
function two_channle_mat_plot(var1::Array{Float64, 2}, var2::Array{Float64, 2})

  seg_height = 0.1
  seg_width = 0.1
  
  xy_co = ([0.1, 0.1], [0.1, 0.3])

  fig_ob = figure()
  ax = fig_ob[:add_subplot](1, 1, 1)
  ax[:set_aspect]("equal")
  for corn in xy_co
      rect = patch.Rectangle((corn[1], corn[2]), seg_width, seg_height, fill = col_mat)
      ax[:add_artist](rect)
  end
  
end