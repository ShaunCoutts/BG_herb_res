# Spatially explicit metabolic and TS herbicide resistance model to test how 
# the conditions required for both metabolic and target site resistance to 
# develope in a populaiton similtaniously.

# tatkes a lower and upper domain limit and resolution and expantion factor for the seedbank 
# returns a tuple of three other tuples
# 1: above ground values of g being evaluated 
# 2: seedabank values of g being evaluated
# 3: indexes of above ground values of g in the tuple of seedbank values of g
function g_points_builder(low_point = -10, upper_point = 10, res = 0.5, seed_expand = 3)
  above_ground_eval = (low_point : res : upper_point)
  seed_eval = (low_point * seed_expand : res : upper_point * seed_expand) 
  above_ground_ind = tuple([findfirst(seed_eval, x) for x in above_ground_eval]...)
  
  return (above_ground_eval, seed_eval, above_ground_ind)
end

g_points_builder()


# produces a matrix of seed dispersal probabilites from every location evaluated to
# every other location in a 1D landscape. Based on a fitted kenel from Colbach and Sache (2001)
# using inplace modificaiton of an exiting matrix disp_mat
function seed_disp_mat_builder_1D!(disp_mat::Array{Float64, 2}; res = 1, pro_short = 0.48, mean_dist_short = 0.58, 
  pro_seeds_to_mean_short = 0.44, mean_dist_long = 1.65, pro_seeds_to_mean_long = 0.39)
  
    
  shape_short = 1 / (log(1 - pro_seeds_to_mean_short) + 1)
  shape_long = 1 / (log(1 - pro_seeds_to_mean_long) + 1)

  scale_short = (shape_short - 1) / (shape_short * mean_dist_short ^ shape_short) 
  scale_long = (shape_long - 1) / (shape_long * mean_dist_long ^ shape_long) 
 
  all_locs = [0; collect(res * 0.5 : res : floor(size(disp_mat)[1] * 4 * res))] # make an intial landcape larger than needed to normalise the dispersl kernel 
  raw_kernel = zeros(length(all_locs))
  # do the first set of calcualtions for the first location x to every other
  d_raw_0 = 0.0 # the weibul is 0 at x = 0, notice second term in equatioon below 
  dist = all_locs[2]
  d_raw_1 = pro_short * scale_short * dist ^ (shape_short - 2) * 
    exp(-scale_short * dist ^ shape_short) + (1 - pro_short) * scale_long * 
    dist ^ (shape_long - 2) * exp(-scale_long * dist ^ shape_long)
  #trapazoidal intergrate to get total density from dist = 0 to dist = 0.5 * res, assume all seeds go to dist_0
  raw_kernel[1] = min(d_raw_0, d_raw_1) * res * 0.5 + res * 0.25 * abs(d_raw_0 - d_raw_1)
  
  d_raw_0 = d_raw_1
 
  @fastmath for x in 2:(length(all_locs) - 1)
    dist = all_locs[x + 1]
    d_raw_1 = pro_short * scale_short * dist ^ (shape_short - 2) * 
      exp(-scale_short * dist ^ shape_short) + (1 - pro_short) * scale_long * 
      dist ^ (shape_long - 2) * exp(-scale_long * dist ^ shape_long)
    #trapazoidal intergrate to get total density from dist_0 to dist_1, assume all seeds go to dist_0
    raw_kernel[x] = min(d_raw_0, d_raw_1) * res + res * 0.5 * abs(d_raw_0 - d_raw_1)
    
    d_raw_0 = d_raw_1
  end
    
  #find the normalising constant
  norm_const = sum(raw_kernel)
  
  disp_mat[1, :] = raw_kernel[1:size(disp_mat)[1]] / (norm_const * 2)
  
  #use this normalised kernel to fill the rest of the pairwise distances 
  disp_mat[2:end, 1] = disp_mat[1, 2:end]
  for i in 2:size(disp_mat)[1]
    disp_mat[i, i] = disp_mat[1, 1] * 2
    disp_mat[i, (i + 1):end] = disp_mat[1, 2:(end - i + 1)] # fill upper triangle row wise 
    disp_mat[(i + 1):end, i] = disp_mat[1, 2:(end - i + 1)] # fill lower triangle col wise
  end
  disp_mat[end, end] = disp_mat[1, 1]
  
  return nothing
end 

# produces a matrix of seed dispersal probabilites from every location evaluated to
# every other location in a 2D landscape. Based on a fitted kenel from Colbach and Sache (2001)
# using inplace modificaiton of an exiting matrix disp_mat
function seed_disp_mat_builder_2D!(disp_mat::Array{Float64, 2}; res = 1, pro_short = 0.48, mean_dist_short = 0.58, 
  pro_seeds_to_mean_short = 0.44, mean_dist_long = 1.65, pro_seeds_to_mean_long = 0.39)
  
    
  shape_short = 1 / (log(1 - pro_seeds_to_mean_short) + 1)
  shape_long = 1 / (log(1 - pro_seeds_to_mean_long) + 1)

  scale_short = (shape_short - 1) / (shape_short * mean_dist_short ^ shape_short) 
  scale_long = (shape_long - 1) / (shape_long * mean_dist_long ^ shape_long) 
  
  area_correction = 1 / (2 * pi) 
 
  # set up a vector of points to evaluate the kernel on, landscape and is used for the normalisation also a template for all othe dispersal events
  all_locs = collect(res * 0.5 : res : floor(sqrt(size(disp_mat)[1]) * 3 * res)) # make an intial landcape larger than needed to normalise the dispersl kernel 
  raw_kernel = zeros(length(all_locs), length(all_locs)) 
  origin = (convert(Int, ceil(length(all_locs) * 0.5)), convert(Int, ceil(length(all_locs) * 0.5)))
  #first square need all four corners 
  for y in 1:length(all_locs)
    for x in 1:length(all_locs)
      dist_0 = sqrt((all_locs[origin[1]] - all_locs[x] - res * 0.5) ^ 2 + (all_locs[origin[2]] - all_locs[y] - res * 0.5) ^ 2) 
      p_0 = area_correction * (pro_short * scale_short * dist_0 ^ (shape_short - 2) * 
	exp(-scale_short * dist_0 ^ shape_short) + (1 - pro_short) * scale_long * 
	dist_0 ^ (shape_long - 2) * exp(-scale_long * dist_0 ^ shape_long))
      
      dist_1 = sqrt((all_locs[origin[1]] - all_locs[x] + res * 0.5) ^ 2 + (all_locs[origin[2]] - all_locs[y] - res * 0.5) ^ 2) 
      p_1 = area_correction * (pro_short * scale_short * dist_1 ^ (shape_short - 2) * 
	exp(-scale_short * dist_1 ^ shape_short) + (1 - pro_short) * scale_long * 
	dist_1 ^ (shape_long - 2) * exp(-scale_long * dist_1 ^ shape_long))
	
      dist_2 = sqrt((all_locs[origin[1]] - all_locs[x] - res * 0.5) ^ 2 + (all_locs[origin[2]] - all_locs[y] + res * 0.5) ^ 2) 
      p_2 = area_correction * (pro_short * scale_short * dist_2 ^ (shape_short - 2) * 
	exp(-scale_short * dist_2 ^ shape_short) + (1 - pro_short) * scale_long * 
	dist_2 ^ (shape_long - 2) * exp(-scale_long * dist_2 ^ shape_long))
      
      dist_3 = sqrt((all_locs[origin[1]] - all_locs[x] + res * 0.5) ^ 2 + (all_locs[origin[2]] - all_locs[y] + res * 0.5) ^ 2) 
      p_3 = area_correction * (pro_short * scale_short * dist_3 ^ (shape_short - 2) * 
	exp(-scale_short * dist_3 ^ shape_short) + (1 - pro_short) * scale_long * 
	dist_3 ^ (shape_long - 2) * exp(-scale_long * dist_3 ^ shape_long))
      
      raw_kernel[y, x] = 0.25 * res^2 * (p_0 + p_1 + p_2 + p_3) 
    end
  end
    
  raw_kernel /= sum(raw_kernel)
  # use the template to to fill in the from/to dispersal matrix
  edge_offset = convert(Int, sqrt(size(disp_mat)[1]))
  row_count = 1
  for from_y in 1:edge_offset
    for from_x in 1:edge_offset
    
      col_count = 1
      #shift from point to a plane centered on the origin
      cent_y = origin[1] - (from_y + edge_offset)
      cent_x = origin[2] - (from_x + edge_offset)
	    
      for y in 1:edge_offset
	for x in 1:edge_offset
	  y_tar = y + edge_offset
	  x_tar = x + edge_offset
	  #shift target to same relative position on new origin centered plane  
	  y_shifted = y_tar + cent_y
	  x_shifted = x_tar + cent_x
	  disp_mat[row_count, col_count] = raw_kernel[y_shifted, x_shifted]
	  col_count += 1
	end
      end
      
      row_count += 1
    
    end
  end
 
  return nothing
end 

# produces a matrix of pollen dispersal probabilites from every location evaluated to
# every other location in a 1D landscape. Based on a fitted logistic kenel from Klein et al 2006
# using inplace modificaiton of an exiting matrix disp_mat
function pollen_disp_mat_builder_1D!(disp_mat::Array{Float64, 2}; res = 1, a = 32.3, c =3.32)

  all_locs = [0; collect(res * 0.5 : res : floor(size(disp_mat)[1] * 5 * res))] # make an intial landcape larger than needed to normalise the dispersl kernel 
  raw_kernel = zeros(length(all_locs))
  
  const_kernel_term = c / (a ^ (2 / c) * gamma(2 / c) * gamma(1 - (2 / c)))
  # do the first set of calcualtions for the first location x to every other
  d_raw_0 = const_kernel_term #note from the equation below that the second term (that includes distance) = 1 when dist = 0
  dist = all_locs[2]
  d_raw_1 = const_kernel_term * (1 + ((dist ^ c) / a)) ^ (-1)
  #trapazoidal intergrate to get total density from dist_0 to dist_1, assume all seeds go to dist_0
  raw_kernel[1] = min(d_raw_0, d_raw_1) * res * 0.5 + res * 0.25 * abs(d_raw_0 - d_raw_1)
  
  d_raw_0 = d_raw_1

  @fastmath for x in 2:(length(all_locs) - 1)
    dist = all_locs[x + 1]
    d_raw_1 = const_kernel_term * (1 + ((dist ^ c) / a)) ^ (-1)
    #trapazoidal intergrate to get total density from dist_0 to dist_1, assume all seeds go to dist_0
    raw_kernel[x] = min(d_raw_0, d_raw_1) * res + res * 0.5 * abs(d_raw_0 - d_raw_1)
    
    d_raw_0 = d_raw_1
  end
 
  #find the normalising constant
  norm_const = sum(raw_kernel)
  
  disp_mat[1, :] = raw_kernel[1:size(disp_mat)[1]] / (norm_const * 2)
  
  #use this normalised kernel to fill the rest of the pairwise distances 
  disp_mat[2:end, 1] = disp_mat[1, 2:end]
  for i in 2:size(disp_mat)[1]
    disp_mat[i, i] = disp_mat[1, 1] * 2
    disp_mat[i, (i + 1):end] = disp_mat[1, 2:(end - i + 1)] # fill upper triangle row wise 
    disp_mat[(i + 1):end, i] = disp_mat[1, 2:(end - i + 1)] # fill lower triangle col wise
  end
 disp_mat[end, end] = disp_mat[1, 1]
  
  return nothing
end

# produces a matrix of pollen dispersal probabilites from every location evaluated to
# every other location in a 2D landscape. Based on a fitted logistic kenel from Klein et al 2006
# using inplace modificaiton of an exiting matrix disp_mat
function pollen_disp_mat_builder_2D!(disp_mat::Array{Float64, 2}; res = 1, a = 32.3, c =3.32)
  
  const_kernel_term = c / (2 * pi * a ^ (2 / c) * gamma(2 / c) * gamma(1 - (2 / c)))
  # set up a vector of points to evaluate the kernel on, landscape and is used for the normalisation also a template for all othe dispersal events
  all_locs = collect(res * 0.5 : res : floor(sqrt(size(disp_mat)[1]) * 3 * res)) # make an intial landcape larger than needed to normalise the dispersl kernel 
  raw_kernel = zeros(length(all_locs), length(all_locs)) 
  origin = (convert(Int, ceil(length(all_locs) * 0.5)), convert(Int, ceil(length(all_locs) * 0.5)))
  #first square need all four corners 
  for y in 1:length(all_locs)
    for x in 1:length(all_locs)
      dist_0 = sqrt((all_locs[origin[1]] - all_locs[x] - res * 0.5) ^ 2 + (all_locs[origin[2]] - all_locs[y] - res * 0.5) ^ 2) 
      p_0 = const_kernel_term * (1 + ((dist_0 ^ c) / a)) ^ (-1)
      
      dist_1 = sqrt((all_locs[origin[1]] - all_locs[x] + res * 0.5) ^ 2 + (all_locs[origin[2]] - all_locs[y] - res * 0.5) ^ 2) 
      p_1 = const_kernel_term * (1 + ((dist_1 ^ c) / a)) ^ (-1)
	
      dist_2 = sqrt((all_locs[origin[1]] - all_locs[x] - res * 0.5) ^ 2 + (all_locs[origin[2]] - all_locs[y] + res * 0.5) ^ 2) 
      p_2 = const_kernel_term * (1 + ((dist_2 ^ c) / a)) ^ (-1)
      
      dist_3 = sqrt((all_locs[origin[1]] - all_locs[x] + res * 0.5) ^ 2 + (all_locs[origin[2]] - all_locs[y] + res * 0.5) ^ 2) 
      p_3 = const_kernel_term * (1 + ((dist_3 ^ c) / a)) ^ (-1)
      
      raw_kernel[y, x] = 0.25 * res^2 * (p_0 + p_1 + p_2 + p_3) 
    end
  end
    
  raw_kernel /= sum(raw_kernel)
  # use the template to to fill in the from/to dispersal matrix
  edge_offset = convert(Int, sqrt(size(disp_mat)[1]))
  row_count = 1
  for from_y in 1:edge_offset
    for from_x in 1:edge_offset
    
      col_count = 1
      #shift from point to a plane centered on the origin
      cent_y = origin[1] - (from_y + edge_offset)
      cent_x = origin[2] - (from_x + edge_offset)
	    
      for y in 1:edge_offset
	for x in 1:edge_offset
	  y_tar = y + edge_offset
	  x_tar = x + edge_offset
	  #shift target to same relative position on new origin centered plane  
	  y_shifted = y_tar + cent_y
	  x_shifted = x_tar + cent_x
	  disp_mat[row_count, col_count] = raw_kernel[y_shifted, x_shifted]
	  col_count += 1
	end
      end
      
      row_count += 1
    
    end
  end
 
  return nothing
end 


# produces a matrix of seed dispersal probabilites from every location evaluated to
# every other location in a 2D landscape
function seed_disp_mat_builder_2D(size = 100, res = 1)
  all_locs = collect(combinations(collect(0 : res : size), 2)) 
  #make a matrix of size x size to hold the results
  for dist in all_dists
    
  end
end 




function seed_disp_mat_builder_1D(landscape; )


# Takes the current state of the population -> next state of the population  
function single_iter(current_pop::Array{Float64, 2}, next_pop::Array{Float64, 2}; seed_disp_mat::Array{Float64, 2}, 
  pollen_disp_mat::Array{Float64, 2})
 
 
 return next_pop 
end


function multi_iter(num_iter = 10, landscape_size = 10, space_res = 1)

  seed_disp_mat = zeros(convert(Int32, ((landscape_size / space_res) + 1) ^ 2), 
    convert(Int32, ((landscape_size / space_res) + 1) ^ 2))
    
  @time seed_disp_mat_builder_1D!(seed_disp_mat, res = space_res)
  
  pollen_disp_mat = zeros(convert(Int32, (landscape_size / space_res) + 1), 
    convert(Int32, (landscape_size / space_res) + 1))
  
  @time pollen_disp_mat_builder_1D!(pollen_disp_mat, res = space_res)
  
  
  @time seed_disp_mat_builder_2D!(seed_disp_mat, res = space_res, pro_short = 0.48, mean_dist_short = 0.58, 
    pro_seeds_to_mean_short = 0.44, mean_dist_long = 1.65, pro_seeds_to_mean_long = 0.39)
  


end