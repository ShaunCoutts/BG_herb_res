# populatoin model for TSR, desinged to translate actions to a
# population trjectory, which can then be turned into rewards

# make some constants so refrecing to arrays and tuples is clearer
const MATG = 1;
const PATG = 2;
const G1 = 3;
const G2 = 4;

##########################################################################
# set up functions, run to pre-calc stuff

# make the mixing key. This is actually slightly slower than simply
# looping everything, but is safer as the keys can be reused to ensure 
# that idicies in different arrays and contexts match up, which can be 
# very hard to reson about
function make_TSR_mix_index(NG::Int64)

	mat_G = Array{Int64, 1}(NG * NG)
	pat_G = Array{Int64, 1}(NG * NG)

	g = 1
	for g1 in 1:NG
		for g2 in 1:NG

			mat_G[g] = g1
			pat_G[g] = g2

			g += 1

		end
	end

	# also genreate a set of keys for refrencesing G1 = [RR, Rr, rr] 
	# and G2 = [AA, Aa, aa]

	NG1 = convert(Int64, sqrt(NG))
	NG2 = convert(Int64, sqrt(NG))

	G1_ind = repeat(1:NG1, inner = NG2)
	G2_ind = repeat(1:NG2, outer = NG1)

	return (mat_G, pat_G, G1_ind, G2_ind) # note order relates to constants above

end

# make the TSR mixing kernel
function make_TSR_kernel(mix_key::Tuple{Array{Int64, 1}, Array{Int64, 1}, 
	Array{Int64, 1}, Array{Int64, 1}})

	# set up two matricies that mix the two TSRs independetly 
	#M = [1.0 0.0 0.0; # RR x RR
	#      0.5 0.5 0.0; # RR x Rr
	#      0.0 1.0 0.0; # RR x rr
	#      0.5 0.5 0.0; # Rr x RR
	#      0.25 0.5 0.25; # Rr x Rr
	#      0.0 0.5 0.5; # Rr x rr
	#      0.0 1.0 0.0; # rr x RR
	#      0.0 0.5 0.5; # rr x Rr
	#      0.0 0.0 1.0] # rr x rr
	# cnvert to arrays of the structure 
	#       PAT
	# MAT |RR|Rr|rr
	# RR
	# Rr  OFF = RR
	# rr
	#       PAT
	# MAT |RR|Rr|rr
	# RR
	# Rr  OFF = Rr
	# rr
	#       PAT
	# MAT |RR|Rr|rr
	# RR
	# Rr  OFF = rr
	# rr
	#
	# M1 mixing for G1, M2 mixing for G2

	M1 = zeros(3, 3, 3)

	M1[:, :, 1] = [1.0 0.5 0.0;
		       0.5 0.25 0.0;
		       0.0 0.0 0.0]

	M1[:, :, 2] = [0.0 0.5 1.0;
		       0.5 0.5 0.5;
		       1.0 0.5 0.0]

	M1[:, :, 3] = [0.0 0.0 0.0;
		       0.0 0.25 0.5;
		       0.0 0.5 1.0]

	M2 = zeros(3, 3, 3)

	M2[:, :, 1] = [1.0 0.5 0.0;
		       0.5 0.25 0.0;
		       0.0 0.0 0.0]

	M2[:, :, 2] = [0.0 0.5 1.0;
		       0.5 0.5 0.5;
		       1.0 0.5 0.0]

	M2[:, :, 3] = [0.0 0.0 0.0;
		       0.0 0.25 0.5;
		       0.0 0.5 1.0]

	# get the number of crosses and targets
	Ncross = length(mix_key[MATG])
	NG = maximum(mix_key[PATG])

	# build the Ncross*Ncross x Noff*Noff joint mixing kernel for G1 and G2
	joint_kernel = Array{Float64, 2}(NG, Ncross)

	for G_off in 1:NG

		cros = 1

		for mat in 1:length(mix_key[G1])
			for pat in 1:length(mix_key[G2])

				joint_kernel[G_off, cros] =	
					M1[mix_key[G1][mat], mix_key[G1][pat], mix_key[G1][G_off]] * 
					M2[mix_key[G2][mat], mix_key[G2][pat], mix_key[G2][G_off]]

					cros += 1 

			end
		end
	end

	# NOTE STRUCTURE OF OUTPUT
	# offspr. (RR,AA)x(RR,AA) (RR,AA)x(RR,Aa) (RR,AA)x(RR,aa) (RR,AA)x(Rr,AA)      
	# RR AA
	# RR Aa
	# RR aa
	# Rr AA
	# Rr Aa
	# Rr aa                 9 x 9^2  offspring ratios      
	# rr AA
	# rr Aa
	# rr aa 

	return joint_kernel

end

# create an intial population vector assuming both TSRs are independent
# we assume the genotypes are well mixed in the seedbank. We have to do 
# this otherwise the simple answer is to flip the seedbank immediatley
function make_int_pop(N1::Float64, N2::Float64, RR::Float64, Rr::Float64,
	AA::Float64, Aa::Float64, mix_key::Tuple{Array{Int64, 1}, 
	Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}})

	# check the arguments are properly constrained
	if (RR + Rr) > 1.0

		error("intial G1 frequencies must sum to < 1")

	end
	if (AA + Aa) > 1.0

		error("intial G2 frequencies must sum to < 1")

	end

	# set up a couple of vectors to loop through
	set_G1 = [RR, Rr, (1 - RR - Rr)]
	set_G2 = [AA, Aa, (1 - AA - Aa)]

	# make the ordered pairs for the population

	SB1 = N1 * set_G1[mix_key[G1]] .* set_G2[mix_key[G2]]
	SB2 = N2 * set_G1[mix_key[G1]] .* set_G2[mix_key[G2]]

	# builds a vector of the pairs
	# SB1 = RR AA | RR Aa | RR aa | Rr AA | Rr Aa | Rr aa | rr AA | rr Aa | rr aa

	return (SB1, SB2)
	
end

# pre-calc the survival that maps the herb action index and genotypeto 
# a survival, where Resistance is dominant
function make_herb_sur_dom(mix_key::Tuple, s0::Float64, pr_ex::Float64,
	sur_herb::Float64)

	# survivals to reference by mk[G1] and mk[G2]
	Hs = [1.0, 1.0, sur_herb]
	# to make herb 1 and 2 have different dominance
	# Hs2 = [1.0, sur_herb, sur_herb]
	# change rest of function accrdingly

	# a set of vectors to hold the survival of each G under each 
	# herb combination, so 4 * 9 = 36 total survivals.
	sur_h0 = ones(length(mix_key[G1])) * s0
	sur_h1 = ((1 - pr_ex) * s0) + 
		(pr_ex * s0 * Hs[mix_key[G1]])
	sur_h2 = ((1 - pr_ex) * s0) + 
		(pr_ex * s0 * Hs[mix_key[G2]])
	sur_h12 = ((1 - pr_ex) * s0) + 
		(pr_ex * s0 * Hs[mix_key[G1]] .* Hs[mix_key[G2]])

	return (sur_h0, sur_h1, sur_h2, sur_h12)

end
# rececive version
function make_herb_sur_rec(mix_key::Tuple, s0::Float64, pr_ex::Float64,
	sur_herb::Float64)

	# survivals to reference by mk[G1] and mk[G2]
	Hs = [1.0, sur_herb, sur_herb]

	# a set of vectors to hold the survival of each G under each 
	# herb combination, so 4 * 9 = 36 total survivals.
	sur_h0 = ones(length(mix_key[G1])) * s0
	sur_h1 = ((1 - pr_ex) * s0) + 
		(pr_ex * s0 * Hs[mix_key[G1]])
	sur_h2 = ((1 - pr_ex) * s0) + 
		(pr_ex * s0 * Hs[mix_key[G2]])
	sur_h12 = ((1 - pr_ex) * s0) + 
		(pr_ex * s0 * Hs[mix_key[G1]] .* Hs[mix_key[G2]])

	return (sur_h0, sur_h1, sur_h2, sur_h12)

end


###########################################################################
# function that do stuff at iteration

# makes the offspring for a given population repreented by mat and pat
# in general pat should be a normalized version of mat
function TSR_par_cross!(mat::Array{Float64, 1}, pat::Array{Float64, 1},
	parents::Array{Float64, 1}, mix_key::Tuple{Array{Int64, 1}, 
	Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}})

	for g in 1:length(mix_key[MATG])

		parents[g] = mat[mix_key[MATG][g]] * 
			pat[mix_key[PATG][g]]

	end

	return nothing

end

# mix seeds between seed bank levels in response to plowing
function plow_inversion!(seedbank_1::Array{Float64, 2}, 
	seedbank_2::Array{Float64, 2}, t::Int64, plow::Bool, 
	inv_frac::Float64)

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

# simulate the population advancing one time step
function onestep!(SB1::Array{Float64, 2}, SB2::Array{Float64, 2}, 
	ab_pop::Array{Float64, 2}, ab_pop_spot::Array{Float64, 2},
	eff_pop::Array{Float64, 1}, mat_pop::Array{Float64, 1}, 
	pat_pop::Array{Float64, 1}, offspr::Array{Float64, 1},
	mix_kern::Array{Float64, 2}, 
	mix_key::Tuple{Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}},
	herb_sur::Tuple{Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}}, 
	inv_frac::Float64, germ_prob::Float64,
	seed_sur::Float64, sur_crop_alt::Float64, sur_spot::Float64,
	fec_dd::Float64, fec_max::Float64,
	herb_act::Int64, crop_act::Int64, plow::Bool, spot_act::Int64, 
	t::Int64)

	# plowing
	plow_inversion!(SB1, SB2, t, plow, inv_frac)
  
	# germination and seed mortality 
	germ_seed_sur!(SB1, SB2, ab_pop, t, germ_prob, seed_sur) 
  
	# survival
	## turnin in to ifelse staments to cover all cropand spot cases (6)
	if crop_act == CROP_FAL
  
		ab_pop[t, :] .= 0.0
  
	elseif crop_act == CROP_WW
    
		ab_pop[t, :] = ab_pop[t, :] .* herb_sur[herb_act]
    
	else
    
		ab_pop[t, :] = (ab_pop[t, :] .* herb_sur[herb_act]) * 
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
	tot_ab_pop = sum(ab_pop_spot[t, :])
  
	# get effective reporduction, taking into account density effects and 
	# demographic costs of resistance
	if tot_ab_pop == 0.0

		eff_pop[:] .= 0.0

	else

		eff_pop[:] = (ab_pop_spot[t, :] ./ (1 + fec_dd * tot_ab_pop)) 
		
	end
 
	# turn ab_pop to frequency dist for paternal distribution
	if tot_ab_pop == 0.0

		pat_pop[:] .= 0.0

	else

		pat_pop[:] = eff_pop ./ sum(eff_pop)

	end
  
	mat_pop[:] = eff_pop * fec_max
  
	# mix seeds over (g1, g2), add the resulting offspring distribution 
	# to seedbank
	TSR_par_cross!(mat_pop, pat_pop, offspr, mix_key);

	SB1[t, :] += mix_kern * offspr
  
	return nothing

end

function one_run!(SB1::Array{Float64, 2}, SB2::Array{Float64, 2}, 
	ab_pop::Array{Float64, 2}, ab_pop_spot::Array{Float64, 2},
	eff_pop::Array{Float64, 1}, mat_pop::Array{Float64, 1}, 
	pat_pop::Array{Float64, 1}, offspr::Array{Float64, 1},
	mix_kern::Array{Float64, 2}, 
	mix_key::Tuple{Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}, Array{Int64, 1}},
	herb_sur::Tuple{Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}, Array{Float64, 1}}, 
	inv_frac::Float64, germ_prob::Float64,
	seed_sur::Float64, sur_crop_alt::Float64, sur_spot::Float64,
	fec_dd::Float64, fec_max::Float64,
	herb_act_seq::Array{Int64, 1}, crop_act_seq::Array{Int64, 1},
	plow_seq::Array{Bool, 1}, spot_act_seq::Array{Int64, 1},	  
	int_SB1::Array{Float64, 1}, int_SB2::Array{Float64, 1}, T::Int64)
  
	# re-set all the values from the previous run that will be 
	# over-written, just to make sure 
	ab_pop[:, :] .= 0.0
  
	# set the intial seed bank
	SB1[1, :] = deepcopy(int_SB1)
	SB2[1, :] = deepcopy(int_SB2)
  
	for t in 2:(T + 1)
  
		SB1[t, :] = deepcopy(SB1[t - 1, :])
		SB2[t, :] = deepcopy(SB2[t - 1, :])
  
		onestep!(SB1, SB2, ab_pop, ab_pop_spot, eff_pop, mat_pop, 
			pat_pop, offspr, mix_kern, mix_key, herb_sur, inv_frac, 
			germ_prob,seed_sur, sur_crop_alt, sur_spot, fec_dd, 
			fec_max, herb_act_seq[t - 1], crop_act_seq[t - 1], 
			plow_seq[t - 1], spot_act_seq[t - 1], t);

	end
  
	return nothing

end



