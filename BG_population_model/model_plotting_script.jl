# plotting the spatial model to see what is going on and start to understand what is happening 
using Colors
#using Gadfly
using PyPlot
using Cairo
using HDF5 JDL #for reading data objects from the disk

# the model output will be complex and large, having to be accesed from disk, so I define a 
# set of handler functions to extract the data 


