#set of functions for herbicide resitance model covering various processes like genotype production and seed bank dynamics
working_loc = '/home/shauncoutts/Dropbox/projects/MHR_blackgrass/BG_model' 
## QUANT_GEN_OFFSPRING_DISTRIBUTION(N_m, N_f, eval_points, additive_variance, offspring_dist_res) 
## produces a matrix where the rows are the distribution of offspring for each evaluated maternal breeding value based on two parent distributions on the breeeding vlaue (N_m and N_f) and 
##a vector of evaluation points on breeding value (each row in returned matrix is the distribution of offspring breeding values returned by each maternal breeding value evaluated), 
## along with parameter for the variance of breeding value for distribtuion of offspring breeding value (which we assume is normal).
## and a resolution to evluate the conditional offspring distribution at
#be careful this functions relies on N_f, N_m, being evaluated on eval_points before being passed to this function so that the indexes all match up
quant_gen_offspring_distribution <- function(N_m, N_f, eval_points, additive_variance, offspring_dist_res){
  additive_sd = sqrt(additive_variance)
  eval_grid = cbind(rep.int(eval_points, length(eval_points)), rep(eval_points, each = length(eval_points)))#make every combination of evaluation points
  index_comb = cbind(rep.int(1:length(eval_points), length(eval_points)), rep(1:length(eval_points), each = length(eval_points)))#make every combination of index on the 1st and 2nd dimention of eval_grid
  N_fathers = N_f / sum(N_f) #turns the distribution of fathers into a frequency distribution
  lower_eval_cond_offspring = 0.5 * eval_grid[1, 1] + 0.5 * eval_grid[1, 2]
  upper_eval_cond_offspring = 0.5 * eval_grid[dim(eval_grid)[1], 1] + 0.5 * eval_grid[dim(eval_grid)[1], 2]
  cond_offspring_eval_points = seq(lower_eval_cond_offspring - additive_sd * 4, upper_eval_cond_offspring + additive_sd * 4, offspring_dist_res) 
  offspring_3D_kernel = sapply(seq_along(index_comb[,1]), FUN = function(x){
    cond_offspring_dist = dnorm(cond_offspring_eval_points, 0.5 * eval_grid[x, 1] + 0.5 * eval_grid[x, 2], additive_sd) #centers the conditional distribtuion of offspring on the breeding value being assesed so that at the extreams the full distribution is still assesed
    cond_offspring_dist = cond_offspring_dist / sum(cond_offspring_dist) #scales the distribtuion so that it sums to one, the assumption being that all offspring produced by a parental combination have to have a breeding value
    cond_offspring_dist * N_fathers[index_comb[x, 1]]  
  })
  offspring_kernel_dims = dim(offspring_3D_kernel)
  summing_grid = matrix(1:(length(eval_points) * length(eval_points)), ncol = length(eval_points), byrow = TRUE)
  sapply(1:offspring_kernel_dims[1], FUN = function(i) apply(summing_grid, MARGIN = 1, FUN = function(x) sum(offspring_3D_kernel[i, x])))
}
#test of quant_gen_offspring_distribution() with dummy input and a fixed known output so I can see if I break this at some point in the future
setwd(working_loc)
test_obj_name <- load('nonspatial_model_test_answer_key.Rdata') #load the test key
eval(parse(text = nonspace_test_answer_key[[1]]$question))#set parameters for the test run
test1 = quant_gen_offspring_distribution(N_m = N_m, N_f = N_f, eval_points = eval_points, additive_variance = additive_variance, offspring_dist_res = offspring_dist_res) #get the output form the current version of the function
#evaluate this aginst the reference answer
#write the output to a named list stored on disk so that this functions output can be tested aginast it at future date if it is changed.
#nonspace_test_answer_key = list()
#nonspace_test_answer_key[[1]] = list(question = 'eval_points = seq(-10, 10, 1.5)\nN_m = N_f = dnorm(eval_points, 0, 2)\nadditive_variance = 0.5\noffspring_dist_res = 0.1\n',
#  answer = test1)
#setwd(working_loc)
#save(nonspace_test_answer_key, file = 'nonspatial_model_test_answer_key.Rdata')

## SURVIVAL()
## produces a distribution (i.e. vector on eval_points) of survivors after herbicide application from a distribution of indivduals that emerge from the seed bank (N_0). 
## eval_points = vector of g values to evaluate the populaiton over
## herb_rate = amount of herbicide applied
## sur_g_min_herb_0 = survival for a minimally NTSR individual when herbicide applicaiton is 0
## g_s50 = value of g at which survival is half sur_g_min_herb_0 when herbicide application is 0
## alpha_s = shape parameter controls how quickley increasing NTSR decreases maxium survival
## zeta = shape parameter that conrtrols how increasing herbicide dose decreases survival for susceptable individuals
## g2resist = slope parameter that translates g score into experienced resistance
## g_res50 = is the resistance score where survival is half survival_g_max, controls how effecive herbicide is. 
survival <- function(N_0, eval_points, herb_rate, sur_g_min_herb_0, g_s50, alpha_s, zeta, g2resist, g_res50){
  survival_g_max = (sur_g_min_herb_0) / (1 + exp(alpha_s * (eval_points - g_s50))) #survival rate for compleatly compleatly resistant individual when herb_rate = 0, can be though of as a survival cost of NTSR
  survival_g_min = (survival_g_max * exp(-zeta * herb_rate)) / (1 + survival_g_max * (exp(-zeta * herb_rate) - 1)) #survial rate for minimmaly resistant individual under herb_rate
  resiatance = eval_points * g2resist  #resistnce score set by eval_points and a slope parameter 
  density_independent_survival = survival_g_min + ((survival_g_max - survival_g_min) / (1 + exp(abs(resiatance - g_res50))))   
  density_independent_survival
}

eval_points = seq(-10, 10, 1.5)
N_0 = dnorm(eval_points, 0, 2)  
herb_rate = 2 
sur_g_min_herb_0 = 0.8 
g_s50 = 100
alpha_s = 0.5
zeta = 0.5
g2resist = 10 
g_res50 = 3

out = survival(N_0 = N_0, eval_points = eval_points, herb_rate = herb_rate, sur_g_min_herb_0 = sur_g_min_herb_0, g_s50 = g_s50, alpha_s = alpha_s, zeta = zeta, g2resist = g2resist, g_res50 = g_res50)
plot(eval_points,out)
#So this does not work, low resistnce indviduals always have better survival reguardless of herbiside rate  

  
  
  
#NOTE TO SELF add non-heritable variance in resitance in the fecundity function so individuals can be resistant through their life time, basically n(g) should be n(g, z) and resistance should then 
#be a function of g and z, so that survival is actually a distribtuion for each element of eval_points do simple version for now.
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#Check the test results  
test_results <- ifelse(all(test1 == nonspace_test_answer_key[[1]]$answer), 'QUANT_GEN_OFFSPRING_DISTRIBUTION() still fine', 'Something you did broke the function QUANT_GEN_OFFSPRING_DISTRIBUTION()')
 
  
  
  
  
  
  
  
  ##area to test things with throw away code#######################################################################################################################
mat = cbind(g_m = rep.int(1:10, length(1:10)), g_f = rep(1:10, each = length(1:10)))
aperm('dim<-' (mat, list(10, 10, 2)), c(2, 1, 3))
