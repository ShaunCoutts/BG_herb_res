# Script to explore phenotyping data for populations we have data for
# goal is to find populations that are very suseptable to all the herbicides and those that 
# are all resistant, then look at population densties of those populations from helen data 
# or just the google drive

setwd('/home/shauncoutts/Dropbox/projects/MHR_blackgrass/Data')
phen_data = read.csv('DS_RES_ForShaun.csv', header = TRUE, stringsAsFactors = FALSE)
crop_data = read.csv('CropData_Shaun.csv', header = TRUE, stringsAsFactors = FALSE) 

# some data clean up
summary(phen_data)
# I want to have one rwo for each field instead of 3 rows fro each field, use field ID
# First take the unique value of co-variates for each field, then add the resistance to each of the three 
# herbicides for each field

phen_list = tapply(seq_along(phen_data$field_name), INDEX = phen_data$field_name, FUN = function(x){
  x_data = phen_data[x, ]
  data.frame(field_name = unique(x_data$field_name)[1], mean_dens = unique(x_data$mean_density_state)[1], farm_name = unique(x_data$farm_name)[1],
    lat_WGS = unique(x_data$latitude_WGS)[1], lon_WGS = unique(x_data$longitude_WGS)[1], n_ATL = x_data[x_data$herb == 'ATL' , 'plants'], 
    n_FEN = x_data[x_data$herb == 'FEN' , 'plants'], n_CYC = x_data[x_data$herb == 'CYC' , 'plants'], dry_wt_ATL = x_data[x_data$herb == 'ATL' , 'dry_wt'], 
    dry_wt_FEN = x_data[x_data$herb == 'FEN' , 'dry_wt'], dry_wt_CYC = x_data[x_data$herb == 'CYC' , 'dry_wt'], sur_ATL = x_data[x_data$herb == 'ATL' , 'survive'], 
    sur_FEN = x_data[x_data$herb == 'FEN' , 'survive'], sur_CYC = x_data[x_data$herb == 'CYC' , 'survive'])
})

phen_df = do.call('rbind', lapply(phen_list, as.data.frame))
row.names(phen_df) <- NULL

# take out some NA's and any feilds where the number of plants in the phenotyping trial is <10
phen_df = phen_df[!is.na(phen_df$n_ATL) & !is.na(phen_df$n_FEN) & !is.na(phen_df$n_CYC) & phen_df$n_ATL > 10 & phen_df$n_FEN > 10 & phen_df$n_CYC > 10, ]
phen_df$p_sur_ATL = phen_df$sur_ATL / phen_df$n_ATL
phen_df$p_sur_FEN = phen_df$sur_FEN / phen_df$n_FEN
phen_df$p_sur_CYC = phen_df$sur_CYC / phen_df$n_CYC
phen_df$min_sur = sapply(seq_along(phen_df$n_ATL), FUN = function(x) min(c(phen_df$p_sur_ATL[x], phen_df$p_sur_FEN[x], phen_df$p_sur_CYC[x])))
phen_df$min_dry_wt = sapply(seq_along(phen_df$n_ATL), FUN = function(x) min(c(phen_df$dry_wt_ATL[x], phen_df$dry_wt_FEN[x], phen_df$dry_wt_CYC[x])))
#also the density class is bound bewteen 0 and 4, rescale [0, 1] so can logit transform
phen_df$rs_dens = phen_df$mean_dens / 4

#now add the crop rotation information
#take the coloumns I want from crop_data
crop_data2 = crop_data[, c('farmname', 'fieldname', 'h.year', 'crop')]
#only take the years 2009 - 2014
crop_data2 = crop_data2[crop_data2$h.year >= 2009 & crop_data2$h.year <= 2014,]
crop_data2 = crop_data2[crop_data2$fieldname != 'NA.', ]
#flatten to make one row per field 
rot_list <- tapply(seq_along(crop_data2$h.year), INDEX = crop_data2$fieldname, FUN = function(x){
  return(rbind(crop_data2$h.year[x], crop_data2$crop[x]))
})

rot_field_names = names(rot_list)
rot_df = data.frame(field_name = rot_field_names, crop_y_2009 = sapply(rot_list, FUN = function(x) x[2, which(x[1, ] == '2009')]), 
  crop_y_2010 = sapply(rot_list, FUN = function(x) x[2, which(x[1, ] == '2010')]), crop_y_2011 = sapply(rot_list, FUN = function(x) x[2, which(x[1, ] == '2011')]), 
  crop_y_2012 = sapply(rot_list, FUN = function(x) x[2, which(x[1, ] == '2012')]), crop_y_2013 = sapply(rot_list, FUN = function(x) x[2, which(x[1, ] == '2013')]), 
  crop_y_2014 = sapply(rot_list, FUN = function(x) x[2, which(x[1, ] == '2014')]), stringsAsFactors = TRUE) 

#join with the phenotype data
phen_rot_df = merge(phen_df, rot_df, by = 'field_name', sort = FALSE, all.x = FALSE)
#make a few summaries of the crop rotation 
phen_rot_df$last2years = phen_rot_df$crop_y_2014 == 'ww' & phen_rot_df$crop_y_2013 == 'ww'
phen_rot_df$last4year =  rowSums(phen_rot_df[,c('crop_y_2014', 'crop_y_2013', 'crop_y_2012', 'crop_y_2011')] == 'ww') / 4  
  
# look at just the populations that have recently been in wheat 3 out of the last 4 year and final year was in winter wheat
phen_l4ww = phen_rot_df[phen_rot_df$last4year >= 0.75, ]
phen_l4ww = phen_ww[phen_ww$crop_y_2014 == 'ww', ]
  
#also look at only the populations thag spent 2013 and 2014 in winter wheat
phen_l2ww = phen_rot_df[phen_rot_df$last2years, ]
  
plot(phen_df)
#more targeted plots to see look for interesitng relationships in data
par(mfrow = c(2, 2))
plot(phen_df$lat_WGS, phen_df$min_sur)
plot(phen_df$min_sur, phen_df$mean_dens)
plot(phen_df$lat_WGS, phen_df$mean_dens)
plot(phen_df$min_sur, qlogis(phen_df$rs_dens))

#same plots with only the ww dominated field
par(mfrow = c(2, 2))
plot(phen_l4ww$lat_WGS, phen_l4ww$min_sur)
plot(phen_l4ww$min_sur, phen_l4ww$mean_dens)
plot(phen_l4ww$lat_WGS, phen_l4ww$mean_dens)
plot(phen_l4ww$min_sur, qlogis(phen_l4ww$rs_dens))

par(mfrow = c(2, 2))
plot(phen_l2ww$lat_WGS, phen_l2ww$min_sur)
plot(phen_l2ww$min_sur, phen_l2ww$mean_dens)
plot(phen_l2ww$lat_WGS, phen_l2ww$mean_dens)
plot(phen_l2ww$min_sur, qlogis(phen_l2ww$rs_dens))

# I want to look for a humped relationship between density and resistance (proxy_min_sur) in the datam which we would expect 
# if there was either over-yeilding or a demographic cost to herbicide resistance  

#fit a gam to the relationship between min_sur (resistance proxy) and logit(density) 
library(mgcv)
gam_fit = gam(qlogis(rs_dens) ~ s(min_sur), data = phen_df)
plot(phen_df$min_sur, qlogis(phen_df$rs_dens))
lines(phen_df$min_sur, gam_fit$fitted.values)
#looks like a very striaght line so either the relationship does not exist or there is too much noise around the data to see a result 

# next goal is to try and use this data to say what we might see in both very resistant and very suceptable fieds  
# try sub-seeting the data to tak only those populations that are most and least suseptable  
phen_extreams = phen_df[phen_df$min_sur > 0.8 | phen_df$min_sur < 0.2,]
par(mfrow = c(2, 2))
plot(phen_extreams$lat_WGS, phen_extreams$min_sur)
plot(phen_extreams$min_sur, phen_extreams$mean_dens)
plot(phen_extreams$lat_WGS, phen_extreams$mean_dens)
plot(phen_extreams$min_sur, qlogis(phen_extreams$rs_dens))

#TODO use the extreams data to set upper and lower limtis to observed number of plants from model at start (i.e. low resistance) and end of run(high resistance)
# a the start of the run we assume the populations will have very low resistance so the lower limit early on (say years 1 - 3) is 0.
# in the upper limit for suseptable populations we can still see mean density classes up to 2.6.
max_den_suc = max(phen_extreams$mean_dens[phen_extreams$min_sur < 0.2]) 
max_den_res = max(phen_extreams$mean_dens[phen_extreams$min_sur > 0.9]) 
min_den_res = min(phen_extreams$mean_dens[phen_extreams$min_sur > 0.9]) 

# however these density classes are not actual numers of plants, but they are realted. to get to a distrbution of the number of plants a 
# a given mean density class represents we can use a burt-force randomisation 

#first a functiont that finds a random matrix with a given density class
# we assume a 1ha field so 5 * 5 = 25 20m x 20m squares that can be in one of 5 density classes 0 - 4 (inclusive)

# some constants to use in the functions
NUM_GRID_SQ = 25
DEN_CLASSES = 0:4

rand_field_gen <- function(mean_dens, acceptable_error, max_iter){
  
  den_vector = sample(DEN_CLASSES, size = NUM_GRID_SQ, replace = TRUE)
  iter = 1
  error_current = abs(mean_dens - mean(den_vector))
  while(iter <= max_iter & error_current >= acceptable_error){
    iter = iter + 1
    #change 1 value at random and chenck if it makes the error less, if so keep the change 
    den_vect_tent = den_vector
    den_vect_tent[sample.int(NUM_GRID_SQ, 1)] = sample(DEN_CLASSES, 1)
    error_tent = abs(mean_dens - mean(den_vect_tent))
    if(error_tent < error_current){
      den_vector = den_vect_tent
      error_current = error_tent
    }
  }
  
  if(error_current > acceptable_error){
    cat('max_iter reached without finding a good solution\n')
    return(paste0('error level achived: ', error_current))
  }else{
    return(den_vector)
  }
}

field_test = rand_field_gen(2.6, 0.01, 10000)

#function to make these randomly generated fields into distrbutions of numbers of plants 
# we use Queenborough et al 2011, Figure 3 for the number of plants these are in ln((plants / m^2) + 1)
# to turn into number of actual plants is (exp(x) - 1) * 400
  
DEN_CLASS_LIMITS = rbind(lower = c(ab = 0, low = 0, med = 69, high = 525, v_high = 1593), 
  upper = c(ab = 69, low = 2121, med = 2900, high = 7800, v_high = 7634))
  
field_2_num <- function(mean_dens, num_rands = 10, acceptable_error = 0.01, max_iter = 10000){
  return(sapply(1:num_rands, FUN = function(k){
    sum(sapply(rand_field_gen(mean_dens = mean_dens, acceptable_error = acceptable_error, max_iter = max_iter), 
      FUN = function(x) runif(1, DEN_CLASS_LIMITS[1, x + 1], DEN_CLASS_LIMITS[2, x + 1])))
    })
  )
}  

hist(field_2_num(4, num_rands = 10000))

# try and find upper and lower limits for populations under no herb and under herb
# take the top and bottom 20% of observations to see the population ranges at the 
# most and least resistant populations, use the phen_l4ww data, has lots of populations
# and they have all been dominated by winter wheat for the last 4 years
summary(phen_l4ww$mean_dens[phen_l4ww$min_sur <= 0.2])
summary(phen_l4ww$mean_dens[phen_l4ww$min_sur >= 0.8])

# assume max mean_dens for a resistant populations is 4 and min is 0.9372
# For a suceptable population assume the min is 0 and max is 2.5780

# turn these numbers into number of plants # min number of plants expected for resistant populations under winter wheat
# minimum number of plants in a suceptable population, assume mean is 0
min_num_plants_suc = field_2_num(0, num_rands = 10000, acceptable_error = 0.01)
min_suc_95 = quantile(min_num_plants_suc, prob = c(0.025, 0.975))
# minmimum number of plants in a suceptable populaiton	
# 0.025 = 667 and 0.975 = 1,055

#maximum number of plants in suceptable population 
max_num_plants_suc = field_2_num(2.578, num_rands = 10000, acceptable_error = 0.02)
max_suc_95 = quantile(max_num_plants_suc, prob = c(0.025, 0.975))
# max number of plants expected for suceptable population under herbicide
# 0.025 = 58,875 and 0.975 = 89,714

# min number of plants expected for resistant populations under winter wheat
min_num_plants_res = field_2_num(0.9372, num_rands = 10000, acceptable_error = 0.02)
min_res_95 = quantile(min_num_plants_res, prob = c(0.025, 0.975))
#minimum number of plants in a resistant population is 
# 0.025 = 16,348 and 0.975 = 33,335

# maximum numebr of plants in resistant population
max_num_plants_res = field_2_num(4, num_rands = 10000, acceptable_error = 0.02)
max_res_95 = quantile(max_num_plants_res, prob = c(0.025, 0.975))
# maximum number of plants in a resistant populiton is 
# 0.025 = 97,941 and 0.975 = 132,003










