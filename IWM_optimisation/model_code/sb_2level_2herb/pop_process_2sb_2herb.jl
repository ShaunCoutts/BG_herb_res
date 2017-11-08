# Population functions for the 2 herbicide, 2 seedbank level model of resistance
using Distributions

#################################################################################################################
##############################################################################################################
# functionst that pre-calculate and set things up

function offspring_dist_setup(g1_vals::Array{Float64, 1}, 
	g2_vals::Array{Float64, 1}, cov_mat::Array{Float64, 2}, 
	mixing_keys::Tuple{Array{Int32, 1}, Array{Int32, 1}})
  
	# array to hold the results, this is going to be messy, 
	# first 2 dimentions will be maternal values of g1 and g2, 
	# next two will be paternal values of g1 and g2
	# then for every quadruple of g values a 2D offspring distribution 
	# will be created over g1 and g2. NOTE this is very senstive to the 
	# size of g1 and g2, so dg has to be quiet course, gd = 1 means 
	# template is 0.68GB, gd = 0.5 means template is 38GB

	len_g1 = size(g1_vals)[1]
	len_par_cross = size(mixing_keys[1])[1]

	#create the combinations of g1, g2 values to evaluate each MVN over
	g_vals_comb = transpose(hcat(g1_vals, g2_vals))

	template = zeros(len_g1, len_par_cross)

	for g in 1:len_par_cross

	    mu = [(g1_vals[mixing_keys[1][g]] + 
		g1_vals[mixing_keys[2][g]]) / 2, 
		(g2_vals[mixing_keys[1][g]] + g2_vals[mixing_keys[2][g]]) / 2]
	    MVN = MvNormal(mu, cov_mat)
	    template[:, g] = pdf(MVN, g_vals_comb)
    
	end
  
	return template
  
end

# Will need to build a 1D vector version that does what I need, 
# which will mean keeping track of indicies. First step is to create a 
# set of indexing keys that relate the 2D maternal distribution
# to the 2D paternal distribution to the 6D offspring distribution.
function make_index_keys(length_g1::Int64, length_g2::Int64)

  tot_len = length_g1 * length_g2
  
  # make Int32, this is going to be a really long vector and memory 
  # is going to be an issue
  maternal_mixing_key = Array{Int32}(tot_len ^ 2) 
  maternal_mixing_key *= 0
  paternal_mixing_key = Array{Int32}(tot_len ^ 2) 
  paternal_mixing_key *= 0
  
  ind_count = 1
  for gm in 1:tot_len
    for gp in 1:tot_len
	    
      maternal_mixing_key[ind_count] = gm
      paternal_mixing_key[ind_count] = gp
      ind_count += 1
	
    end
  end
 
  return (maternal_mixing_key, paternal_mixing_key)
 
end 

# make the reduction in fec with increasing resistance 
function fec_cost_maker(fr::Float64, f0::Float64, g1::Array{Float64, 1}, 
	g2::Array{Float64, 1})

	return 1 ./ (1 + exp(-(f0 - fr * abs(g1) - fr * abs(g2))))

end

# make the survival 
function survial_herb_setup(g1_vals::Array{Float64, 1}, 
	g2_vals::Array{Float64, 1}, pro_exposed::Float64, sur_base::Float64, 
	effect_herb1::Float64, effect_herb2::Float64, prot_g1_herb1::Float64, 
	prot_g2_herb2::Float64)

	# no herbicide
	sur_none = 1 ./ (1 + exp(-(ones(size(g1_vals)[1]) * sur_base)))

	# herbicide 1
	sur_1 = ((1 - pro_exposed) ./ (1 + exp(-(ones(size(g1_vals)[1]) * 
		sur_base)))) + (pro_exposed ./ (1 + exp(-(sur_base - 
		(effect_herb1 - min(effect_herb1, prot_g1_herb1 * g1_vals))))))
  
	# herbicide 2
	sur_2 = ((1 - pro_exposed) ./ (1 + exp(-(ones(size(g1_vals)[1]) * 
		sur_base)))) + (pro_exposed ./ (1 + exp(-(sur_base - 
		(effect_herb2 - min(effect_herb2, prot_g2_herb2 * g2_vals))))))
  
	# both herbicides
	sur_both = ((1 - pro_exposed) ./ (1 + exp(-(ones(size(g1_vals)[1]) * 
		sur_base)))) + pro_exposed * ((1 ./ (1 + exp(-(sur_base - 
		(effect_herb1 - min(effect_herb1, prot_g1_herb1 * g1_vals)))))) .* 
		(1 ./(1 + exp(-(sur_base - (effect_herb2 - min(effect_herb2, 
		prot_g2_herb2 * g2_vals)))))))
    
	return (sur_none, sur_1, sur_2, sur_both)

end
#################################################################################################################
# Functions that do things during timestep

# use the template and a vector of [n_m * phi * n_p], every combination of 
# n_m:n_p to produce the final distribution of seeds over (g1, g2)
function seed_cross(paternal::Array{Float64, 1}, new_seeds::Array{Float64, 1}, 
	dg::Float64, par_mix::Array{Float64, 1}, 
	mixing_keys::Tuple{Array{Int32, 1}, Array{Int32, 1}}, 
	mixing_kernel::Array{Float64, 2})
  
	# construct the vectors for the convolutions 
	par_mix[:] = paternal[mixing_keys[2]] .* new_seeds[mixing_keys[1]]
  
	#do the convolution
	return (mixing_kernel * par_mix) * dg * dg * dg * dg  
  
end

# mix seeds between sssd bank levels in response to plowing
function plow_inversion!(seedbank_1::Array{Float64, 2}, 
	seedbank_2::Array{Float64, 2}, t::Int64, plow::Bool, inv_frac::Float64)

	  if plow
    
	# seeds on the move
	L1_2_L2 = seedbank_1[t, :] * inv_frac
	L2_2_L1 = seedbank_2[t, :] * inv_frac
    
	seedbank_1[t, :] = seedbank_1[t, :] .- L1_2_L2 .+ L2_2_L1
	seedbank_2[t, :] = seedbank_2[t, :] .- L2_2_L1 .+ L1_2_L2
    
	end
  
	return nothing

end

# germination and seed mortality, updates seedbank and creates above ground population
function germ_seed_sur!(seedbank_1::Array{Float64, 2}, 
	seedbank_2::Array{Float64, 2}, ab_pop::Array{Float64, 2}, t::Int64, 
	germ_prob::Float64, seed_sur::Float64)
  
	  seedbank_1[t, :] = seedbank_1[t, :] * seed_sur
	  seedbank_2[t, :] = seedbank_2[t, :] * seed_sur
  
	  ab_pop[t, :] = seedbank_1[t, :] * germ_prob
  
	  seedbank_1[t, :] = seedbank_1[t, :] * (1 - germ_prob)
  
	  return nothing

end


# takes the seed bank and iterates it through one time step to get the 
# seed bank at the start of the next time step, also record above ground 
# population to use for reward function calculation 
function one_step!(seedbank_1::Array{Float64, 2}, 
	seedbank_2::Array{Float64, 2}, ab_pop::Array{Float64, 2}, 
	ab_pop_spot::Array{Float64, 2}, mat_pop::Array{Float64, 1}, 
	pat_pop::Array{Float64, 1}, eff_pop::Array{Float64, 1}, 
	par_mix::Array{Float64, 1}, mixing_keys::Tuple{Array{Int32, 1}, 
	Array{Int32, 1}}, mixing_kernel::Array{Float64, 2}, 
	crop_sur_tup::Tuple{Float64, Float64, Float64}, 
	g_eff_dem::Array{Float64, 1}, herb_sur_tup::Tuple{Array{Float64, 1}, 
	Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}}, 
	crop_act::Int64, spot_act::Int64, plow::Bool, herb_act::Int64, 
	sur_crop_alt::Float64, inv_frac::Float64, germ_prob::Float64, 
	seed_sur::Float64, fec_max::Float64, fec_dd::Float64, 
	sur_spot::Float64, dg::Float64, t::Int64)
    
	# plowing
	plow_inversion!(seedbank_1, seedbank_2, t, plow, inv_frac)
  
	# germination and seed mortality 
	germ_seed_sur!(seedbank_1, seedbank_2, ab_pop, t, germ_prob, seed_sur) 
  
	# survival
	## turnin in to ifelse staments to cover all cropand spot cases (6)
	if crop_act == CROP_FAL
  
		ab_pop[t, :] .= 0.0
  
	elseif crop_act == CROP_WW
    
		ab_pop[t, :] = ab_pop[t, :] .* herb_sur_tup[herb_act]
    
	else
    
		ab_pop[t, :] = (ab_pop[t, :] .* herb_sur_tup[herb_act]) * 
			sur_crop_alt
  
	end
  
	# apply spot control, but store the pre spot control population to 
	# calculate the cost of spot control
	if spot_act == 1

		ab_pop_spot[t, :] = ab_pop[t, :] * sur_spot

	else

		ab_pop_spot[t, :] = ab_pop[t, :] * 1.0

	end

	# post survial above ground pop number
	tot_ab_pop = sum(ab_pop_spot[t, :]) * dg * dg
  
	# get effective reporduction, taking into account density effects and 
	# demographic costs of resistance
	if tot_ab_pop == 0.0

		eff_pop[:] .= 0.0

	else

		eff_pop[:] = (ab_pop_spot[t, :] .* g_eff_dem) * 
			(1.0 / (1 + fec_dd * tot_ab_pop)) 
	end
 
	# turn ab_pop to frequency dist for paternal distribution
	if tot_ab_pop == 0.0

		pat_pop[:] .= 0.0

	else

		pat_pop[:] = eff_pop / (sum(eff_pop) * dg * dg)

	end
  
	mat_pop[:] = eff_pop * fec_max
  
	# mix seeds over (g1, g2), add the resulting offspring distribution 
	# to seedbank
	seedbank_1[t, :] += seed_cross(pat_pop, mat_pop, dg, par_mix, 
		mixing_keys, mixing_kernel) 
  
	return nothing

end

function one_run!(seedbank_1::Array{Float64, 2}, seedbank_2::Array{Float64, 2},
	ab_pop::Array{Float64, 2}, ab_pop_spot::Array{Float64, 2}, 
	mat_pop::Array{Float64, 1}, pat_pop::Array{Float64, 1}, 
	eff_pop::Array{Float64, 1}, par_mix::Array{Float64, 1}, 
	mixing_keys::Tuple{Array{Int32, 1}, Array{Int32, 1}}, 
	mixing_kernel::Array{Float64, 2}, g_eff_dem::Array{Float64, 1}, 
	crop_sur_tup::Tuple{Float64, Float64, Float64}, 
	herb_sur_tup::Tuple{Array{Float64, 1}, Array{Float64, 1}, 
	Array{Float64, 1}, Array{Float64, 1}}, crop_act_seq::Array{Int64, 1}, 
	spot_act_seq::Array{Int64, 1}, plow_seq::Array{Bool, 1}, 
	herb_act_seq::Array{Int64, 1}, time_hor::Int64, 
	int_sbL1::Array{Float64, 1}, int_sbL2::Array{Float64, 1}, 
	sur_crop_alt::Float64, inv_frac::Float64, germ_prob::Float64, 
	seed_sur::Float64, fec_max::Float64, fec_dd::Float64, 
	sur_spot::Float64, dg::Float64)
  
	# re-set all the values from the previous run that will be 
	# over-written, just to make sure 
	seedbank_1[:, :] .= 0.0
	seedbank_2[:, :] .= 0.0
	ab_pop[:, :] .= 0.0
  
	# set the intial seed bank
	seedbank_1[1, :] = deepcopy(int_sbL1)
	seedbank_2[1, :] = deepcopy(int_sbL2)
  
	for i in 2:(time_hor + 1)
  
		seedbank_1[i, :] = deepcopy(seedbank_1[i - 1, :])
		seedbank_2[i, :] = deepcopy(seedbank_2[i - 1, :])
  
		one_step!(seedbank_1, seedbank_2, ab_pop, ab_pop_spot, mat_pop,
			pat_pop, eff_pop, par_mix, mixing_keys, mixing_kernel,
			crop_sur_tup, g_eff_dem, herb_sur_tup, 
			crop_act_seq[i - 1], spot_act_seq[i - 1], 
			plow_seq[i - 1], herb_act_seq[i - 1], sur_crop_alt, 
			inv_frac, germ_prob, seed_sur, fec_max, fec_dd, 
			sur_spot, dg, i)
	    
	end
  
	return nothing

end

