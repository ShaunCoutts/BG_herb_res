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

###########################################################################

function TSR_par_cross(mat::Array{Float64, 1}, pat::Array{Float64, 1},
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
		L1_2_L2 = seedbank_1[:, t] * inv_frac
		L2_2_L1 = seedbank_2[:, t] * inv_frac
    
		seedbank_1[:, t] = seedbank_1[:, t] .- L1_2_L2 .+ L2_2_L1
		seedbank_2[:, t] = seedbank_2[:, t] .- L2_2_L1 .+ L1_2_L2
    
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
function onestep!()


end

