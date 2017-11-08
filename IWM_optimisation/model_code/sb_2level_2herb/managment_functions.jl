# managment and optimisation functions 

# define action indexes that are const through whole script 
const HERB0 = 1;
const HERB1 = 2;
const HERB2 = 3;
const HERB12 = 4;

const CROP_WW = 1;
const CROP_ALT = 2;
const CROP_FAL = 3;

const PLOW0 = 1;
const PLOW = 2;

const SPOT0 = 0;
const SPOT = 1;

const ACT_HERB = 1;
const ACT_CROP = 2;
const ACT_PLOW = 3;
const ACT_SPOT = 4;

const plow_subact = [false, true];

# create the actions space, really just hard coded 
function make_action_space()

	return (
	[HERB0, CROP_WW, PLOW0, SPOT0], 
	[HERB0, CROP_WW, PLOW, SPOT0], 
	[HERB0, CROP_ALT, PLOW0, SPOT0], 
	[HERB0, CROP_ALT, PLOW, SPOT0], 
	[HERB0, CROP_FAL, PLOW0, SPOT0], 
	[HERB0, CROP_FAL, PLOW, SPOT0], 
	[HERB1, CROP_WW, PLOW0, SPOT0], 
	[HERB1, CROP_WW, PLOW, SPOT0], 
	[HERB1, CROP_ALT, PLOW0, SPOT0], 
	[HERB1, CROP_ALT, PLOW, SPOT0], 
	[HERB2, CROP_WW, PLOW0, SPOT0], 
	[HERB2, CROP_WW, PLOW, SPOT0],
	[HERB2, CROP_ALT, PLOW0, SPOT0], 
	[HERB2, CROP_ALT, PLOW, SPOT0], 
	[HERB12, CROP_WW, PLOW0, SPOT0], 
	[HERB12, CROP_WW, PLOW, SPOT0],  
	[HERB12, CROP_ALT, PLOW0, SPOT0],
	[HERB12, CROP_ALT, PLOW, SPOT0],
	[HERB0, CROP_WW, PLOW0, SPOT], 
	[HERB0, CROP_WW, PLOW, SPOT], 
	[HERB0, CROP_ALT, PLOW0, SPOT], 
	[HERB0, CROP_ALT, PLOW, SPOT], 
	[HERB1, CROP_WW, PLOW0, SPOT], 
	[HERB1, CROP_WW, PLOW, SPOT], 
	[HERB1, CROP_ALT, PLOW0, SPOT], 
	[HERB1, CROP_ALT, PLOW, SPOT], 
	[HERB2, CROP_WW, PLOW0, SPOT], 
	[HERB2, CROP_WW, PLOW, SPOT], 
	[HERB2, CROP_ALT, PLOW0, SPOT], 
	[HERB2, CROP_ALT, PLOW, SPOT], 
	[HERB12, CROP_WW, PLOW0, SPOT], 
	[HERB12, CROP_WW, PLOW, SPOT],  
	[HERB12, CROP_ALT, PLOW0, SPOT],
	[HERB12, CROP_ALT, PLOW, SPOT])

end

# create cost space
function make_cost_space(cost_herb_one::Float64, cost_WW::Float64, 
	cost_ALT::Float64, cost_FAL::Float64, cost_plow::Float64)

	herb_cost = zeros(HERB12)
	herb_cost[HERB0] = 0.0
	herb_cost[HERB1] = cost_herb_one
	herb_cost[HERB2] = cost_herb_one
	herb_cost[HERB12] = 2 * cost_herb_one

	crop_cost = zeros(CROP_FAL)
	crop_cost[CROP_WW] = cost_WW
	crop_cost[CROP_ALT] = cost_ALT
	crop_cost[CROP_FAL] = cost_FAL

	plow_cost = zeros(PLOW)
	plow_cost[PLOW0] = 0.0
	plow_cost[PLOW] = cost_plow

	return (herb_cost, crop_cost, plow_cost)

end
# convert a target survival to a value of g, used to get intial condtions
function sur_2_g(targ_sur::Float64, herb_ef::Float64, s0::Float64, 
		 g_pro::Float64)

	return -((log((1 / targ_sur) - 1) - herb_ef + s0) / g_pro)

end

# pre-calc the a vecotor of discount rates to apply at each time step
function make_dis_rate(T::Int64, dis_rate::Float64)
	
	return discounts = dis_rate .^ collect(0:(T - 1)) 

end


# generate a population random seqence of actions
function rand_act_pop(num_acts::Int64, time_hor::Int64, pop_num::Int64)

  return rand(1:num_acts, (pop_num, time_hor))

end

# takes an action seqence and an action space A and turns it into a 
# seqence of [herb_seq, crop_seq, plow_seq].
function act_seq_2_sub_act(A::Tuple, act_seq::Array{Int64, 1})
  
	a = Array{Int64, 2}(size(act_seq)[1], 4)
 
	for i in 1:size(act_seq)[1]
    
		a[i, ACT_HERB] = A[act_seq[i]][ACT_HERB]
		a[i, ACT_CROP] = A[act_seq[i]][ACT_CROP]
		a[i, ACT_PLOW] = A[act_seq[i]][ACT_PLOW]
		a[i, ACT_SPOT] = A[act_seq[i]][ACT_SPOT]

	end

	return a

end

# generate a couple of special case action sequences for testing 
function act_seq_herb0(A::Tuple, time_hor::Int64)

  a = Array{Int64, 2}(time_hor, 4)
  
  a[:, ACT_HERB] .= HERB0
  a[:, ACT_CROP] .= CROP_WW
  a[:, ACT_PLOW] .= PLOW0
  a[:, ACT_SPOT] .= SPOT0
  
  return a

end

function act_seq_herb1(time_hor::Int64)

  a = Array{Int64, 2}(time_hor, 4)
  
  a[:, ACT_HERB] .= HERB1
  a[:, ACT_CROP] .= CROP_WW
  a[:, ACT_PLOW] .= PLOW0
  a[:, ACT_SPOT] .= SPOT0
  
  return a

end

function act_seq_herb2(time_hor::Int64)

  a = Array{Int64, 2}(time_hor, 4)
  
  a[:, ACT_HERB] .= HERB2
  a[:, ACT_CROP] .= CROP_WW
  a[:, ACT_PLOW] .= PLOW0
  a[:, ACT_SPOT] .= SPOT0
  
  return a

end

function act_seq_herb12(time_hor::Int64)

  a = Array{Int64, 2}(time_hor, 4)
  
  a[:, ACT_HERB] .= HERB12
  a[:, ACT_CROP] .= CROP_WW
  a[:, ACT_PLOW] .= PLOW0
  a[:, ACT_SPOT] .= SPOT0
  
  return a

end

function act_seq_plow(time_hor::Int64)

  a = Array{Int64, 2}(time_hor, 4)
  
  a[:, ACT_HERB] .= HERB0
  a[:, ACT_CROP] .= CROP_WW
  a[:, ACT_PLOW] .= PLOW
  a[:, ACT_SPOT] .= SPOT0
  
  return a

end

# reduction for repeated crop
function repeat_penelty(crop_t, crop_t1, penelty)

  if(crop_t == crop_t1)
    
    return(penelty)
    
  else
    
    return 1.0
    
  end

end

# return a vector of economic rerwards after taking a vector of population sizes and actions 
function economic_reward(N::Array{Float64, 1}, crop::Array{Int64, 1}, 
	Y0::Float64, slope::Float64, alt_Y::Float64, rep_penelty::Float64)

	econ_reward = zeros(length(N))

	# yeild in first year
	if crop[1] == CROP_WW

		econ_reward[1] = (Y0 - slope * N[1]) 

	elseif crop[1] == CROP_ALT

		econ_reward[1] = alt_Y

	elseif crop[1] == CROP_FAL

		econ_reward[1] = 0.0

	else

		econ_reward[i] = NaN

	end

	# following years
	for i in 2:length(econ_reward)

		if crop[i] == CROP_WW

			econ_reward[i] = (Y0 - slope * N[i]) *  
				repeat_penelty(crop[i - 1], crop[i], rep_penelty)

		elseif crop[i] == CROP_ALT

			econ_reward[i] = alt_Y * repeat_penelty(crop[i - 1], 
				crop[i], rep_penelty)

		elseif crop[i] == CROP_FAL

			econ_reward[i] = 0.0

		else

			econ_reward[i] = NaN 
		end
	end

	return econ_reward

end

# return a vector of social rewards after taking vector of populations 
function social_reward(N::Array{Float64, 1}, scale::Float64, shape::Float64)

  return scale * exp(-shape * N)

end

# return a vector of costs taking a matirx of actions (rows = t, colunm = acttions)
function costs(actions::Array{Int64, 2}, cost_space::Tuple{Array{Float64, 1}, 
	Array{Float64, 1}, Array{Float64, 1}}, cost_spot::Float64,
        N_pre_spot::Array{Float64, 1})

	T = size(actions)[1]
	cost_vect = zeros(T)
	spot_cost = 0.0
	for t in 1:T

		if actions[t, ACT_SPOT] == 0

			spot_cost = 0.0

		else

			spot_cost = cost_spot * N_pre_spot[t]

		end

		cost_vect[t] = cost_space[ACT_HERB][actions[t, ACT_HERB]] + 
			cost_space[ACT_CROP][actions[t, ACT_CROP]] + 
			cost_space[ACT_PLOW][actions[t, ACT_PLOW]] +
			spot_cost

	end

	return cost_vect

end

# calcualte the reward from a single run
function reward_total(ab_pop::Array{Float64, 2}, 
	ab_pop_spot::Array{Float64, 2}, dis_rates::Array{Float64, 1}, 
	Y0::Float64, slope::Float64, alt_Y::Float64, rep_penelty::Float64, 
	act_seq::Array{Int64, 2}, cost_space::Tuple{Array{Float64, 1}, 
	Array{Float64, 1}, Array{Float64, 1}}, cost_spot::Float64,
	dg::Float64)
  
	tot_ab_pop = vcat(sum(ab_pop, 2)...) * dg * dg
	tot_ab_pop = tot_ab_pop[2:end]

	tot_ab_pop_spot = vcat(sum(ab_pop_spot, 2)...) * dg * dg
	tot_ab_pop_spot = tot_ab_pop_spot[2:end]

	raw_reward = economic_reward(tot_ab_pop_spot, act_seq[:, ACT_CROP],
			Y0, slope, alt_Y, rep_penelty) - 
		costs(act_seq, cost_space, cost_spot, tot_ab_pop) 

	dis_reward = dis_rates .* raw_reward

	return sum(dis_reward)

end

# evaluate a population of action sequences
function eval_act_seq_pop(pop_act_seq::Array{Int64, 2}, 
	act_space::Tuple, seedbank_1::Array{Float64, 2}, 
	seedbank_2::Array{Float64, 2}, ab_pop::Array{Float64, 2}, 
	ab_pop_spot::Array{Float64, 2}, mat_pop::Array{Float64, 1}, 
	pat_pop::Array{Float64, 1}, eff_pop::Array{Float64, 1}, 
	par_mix::Array{Float64, 1}, 
	mix_keys::Tuple{Array{Int32, 1}, Array{Int32, 1}}, 
  	mix_kernel::Array{Float64, 2}, fit_cost::Array{Float64, 1}, 
	crop_sur_tup::Tuple{Float64, Float64, Float64}, 
  	herb_sur_tup::Tuple{Array{Float64, 1}, Array{Float64, 1}, 
	Array{Float64, 1}, Array{Float64, 1}}, T::Int64, 
	int_sb1::Array{Float64, 1}, int_sb2::Array{Float64, 1}, 
	sur_crop_alt::Float64, inv_frac::Float64, germ_prob::Float64, 
	seed_sur::Float64, fec_max::Float64, fec_dd::Float64, 
	sur_spot::Float64, dg::Float64, dis_rates::Array{Float64, 1}, 
	Y0::Float64, Y_slope::Float64, Y_ALT::Float64, 
	cost_space::Tuple{Array{Float64, 1}, Array{Float64, 1}, 
	Array{Float64, 1}}, rep_pen::Float64, cost_spot::Float64)

	act_pop_size = size(pop_act_seq)[1]
	rewards = zeros(act_pop_size)
    
	for i in 1:act_pop_size

		# get action sequences
		sub_acts = act_seq_2_sub_act(act_space, pop_act_seq[i, :])
		crop_act_seq = sub_acts[:, ACT_CROP]
		spot_act_seq = sub_acts[:, ACT_SPOT]
		plow_seq = plow_subact[sub_acts[:, ACT_PLOW]]
		herb_act_seq = sub_acts[:, ACT_HERB]

		# run under act_seq_pop[i, :] to get Ns
		one_run!(seedbank_1, seedbank_2, ab_pop, ab_pop_spot, mat_pop,
			 pat_pop, eff_pop, par_mix, mix_keys, mix_kernel,
			 fit_cost, crop_sur_tup, herb_sur_tup, crop_act_seq, 
			 spot_act_seq, plow_seq, herb_act_seq, T, 
			 int_sb1, int_sb2, sur_crop_alt, inv_frac, germ_prob,
			 seed_sur, fec_max, fec_dd, sur_spot, dg)

		# get the rewards given act_seq_pop[i, :]
		rewards[i] = reward_total(ab_pop, ab_pop_spot, dis_rates, 
			Y0, Y_slope, Y_ALT, rep_pen, sub_acts, cost_space, 
			cost_spot, dg)

	end
 
  return rewards

end

# tournament selection to get indicies of the survivours 
function tourn_select(rewards::Array{Float64, 1})

	# set up the pairs at random 
	N = length(rewards)
	rand_inds = sample(1:N, N, replace = false)

	pop1 = rand_inds[1:convert(Int64, floor(N / 2))]
	pop2 = rand_inds[(convert(Int64, floor(N / 2)) + 1):N]

	# find the winners
	winners = ifelse(rewards[pop1] .> rewards[pop2], pop1, pop2)

	return winners

end

# do the crosses and mutation
function cross_mut(pop_1::Array{Int64, 1}, pop_2::Array{Int64, 1},
	mut::Float64, num_acts::Int64)

	# find break point
	swap_ind = rand(2:(length(pop_1) - 1))

	# make the swap
	child = vcat(pop_1[1:swap_ind], pop_2[(swap_ind + 1):end])

	# do the mutation and return
	return ifelse(rand(length(child)) .> mut, child, 
		rand(1:num_acts, length(child)))

end


# next generation of action sequences
# note winner is of length size(act_pop_t)[1] / 2
# and size(act_pop_t)[1] is even
function next_gen!(act_pop_t::Array{Int64, 2}, act_pop_t1::Array{Int64, 2},
	winners::Array{Int64, 1}, mut::Float64, num_acts::Int64)
	
	W = length(winners)
	# put the parents in the next 
	for i in 1:W
	
		act_pop_t1[i, :] = act_pop_t[winners[i], :]

	end
	
	# set up the pairs to make babies
	p1 = winners
	p2 = vcat(winners[2:end], winners[1])

	# make the crosses
	for i in 1:W

		act_pop_t1[W + i, :] = cross_mut(act_pop_t[p1[i], :], 
			act_pop_t[p2[i], :], mut, num_acts)

	end

	return nothing

end

# One GA run to find good managment sequences
function GA_solve(T::Int64, pop_size::Int64, num_gen::Int64, 
	cost_herb_one::Float64, cost_WW::Float64, cost_ALT::Float64, 
	cost_FAL::Float64, cost_plow::Float64, cost_spot::Float64, 
	sur_crop_alt::Float64, low_g::Float64, up_g::Float64, dg::Float64, 
	off_sd::Float64, off_cv::Float64, int_N::Float64, int_sd::Float64, 
	int_cv::Float64, int_g1::Float64, int_g2::Float64, inv_frac::Float64, 
	germ_prob::Float64, seed_sur::Float64, fec_max::Float64, 
	fec_dd::Float64, sur_spot::Float64, dis_rate::Float64, Y0::Float64,
	Y_slope::Float64,Y_ALT::Float64, pro_exposed::Float64, 
	sur_base::Float64, rep_pen::Float64, effect_herb1::Float64, 
	effect_herb2::Float64, prot_g1_herb1::Float64, prot_g2_herb2::Float64, 
	fr::Float64, f0::Float64, mut::Float64)

	# hold the population actions 
	pop_list = Array{Array, 1}(num_gen)
	for g in 1:num_gen

		pop_list[g] = Array{Int64, 2}(pop_size, T)

	end
	
	# hold the rewards for each action sequence 
	reward_list = Array{Array, 1}(num_gen)
	for g in 1:num_gen

		reward_list[g] = zeros(pop_size)

	end

	# create a set of spaces and holding arrays
	A = make_action_space()
	C = make_cost_space(cost_herb_one, cost_WW, cost_ALT, cost_FAL,
		cost_plow)

	# make the intial g_val vectors note the repeated values to 
	# encode a 2D vector in 1D
	g_vals = collect(low_g : dg : up_g)
	len_g = size(g_vals)[1]
	g1_vals = repeat(g_vals, inner = len_g)
	g2_vals = repeat(g_vals, outer = len_g)

	Ng = length(g1_vals)

	# pre-calc the effet of managment actions 
	crop_sur_tup = (1.0, sur_crop_alt, 0.0)
	spot_sur_tup = (0.0, 1.0)

	herb_sur_tup = survial_herb_setup(g1_vals, g2_vals, pro_exposed, 
		sur_base, effect_herb1, effect_herb2, prot_g1_herb1, 
		prot_g2_herb2)

	fit_cost = fec_cost_maker(fr, f0, g1_vals, g2_vals)

	dis_rates = make_dis_rate(T, dis_rate)

	# mixing kernels and indexing keys
	mix_keys = make_index_keys(len_g, len_g)
	# define resistance trait values and their co-var
	cov_mat = [off_sd off_cv;
		   off_cv off_sd]
	mix_kernel = offspring_dist_setup(g1_vals, g2_vals, cov_mat, 
		mix_keys)

	# seedbanks T by Ng in size
	SB1 = zeros(T + 1, Ng) 
	SB2 = zeros(T + 1, Ng) 
	# above ground populations 
	ab_pop = zeros(T + 1, Ng)
	ab_pop_spot = zeros(T + 1, Ng)
	# temp holding, rewritten eacbh time step
	mat_pop = zeros(Ng)
	pat_pop = zeros(Ng)
	eff_pop = zeros(Ng)
	par_mix = zeros(length(mix_keys[1]))

	# intial populations 
	int_cov = [int_sd int_cv; 
		   int_cv int_sd]
	int_mu = [int_g1; int_g2]
	int_dist = MvNormal(int_mu, int_cov);
	int_sb1 = pdf(int_dist, transpose(hcat(g1_vals, g2_vals))) * int_N  
	int_sb2 = deepcopy(int_sb1) 

	# get number of acts, it is used a few times
	N_acts = length(A)

	# make the first populaiton random sequences
	pop_list[1][:, :] = rand_act_pop(N_acts, T, pop_size)
	
	# evaluate the sequences 
	for g in 1:(num_gen - 1)
		
		# evaluate each sequence
		reward_list[g][:] = eval_act_seq_pop(pop_list[g], A,
			SB1, SB2, ab_pop, ab_pop_spot, mat_pop,
			pat_pop, eff_pop, par_mix, mix_keys, mix_kernel,
			fit_cost, crop_sur_tup, herb_sur_tup, T, int_sb1, 
			int_sb2, sur_crop_alt, inv_frac, germ_prob, 
			seed_sur, fec_max, fec_dd, sur_spot, dg, dis_rates,
			Y0, Y_slope, Y_ALT, C, rep_pen, cost_spot)

		# indicies of sequences that performed well  
		win_ind = tourn_select(reward_list[g])
		
		# make the next generation of action sequences through 
		# cross and mutation
		next_gen!(pop_list[g], pop_list[g + 1], win_ind, mut, N_acts)

	end

	# evaluate final sequence
	reward_list[end][:] = eval_act_seq_pop(pop_list[end], A,
		SB1, SB2, ab_pop, ab_pop_spot, mat_pop,
		pat_pop, eff_pop, par_mix, mix_keys, mix_kernel,
		fit_cost, crop_sur_tup, herb_sur_tup, T, int_sb1, 
		int_sb2,sur_crop_alt, inv_frac, germ_prob, 
		seed_sur, fec_max, fec_dd, sur_spot, dg, dis_rates, Y0, 
		Y_slope, Y_ALT, C, rep_pen, cost_spot)
	
	# put the parameter values in a Dict for easy use later in 
	# plotting
	par_dict = Dict(:int_g1 => int_g1, 
			:int_g2 => int_g2,
			:int_N => int_N, 
			:int_sd => int_sd, 
			:int_cv => int_cv, 
			:p_ex => pro_exposed, 
			:s0 => sur_base,
			:eff_h1 => effect_herb1, 
			:eff_h2 => effect_herb2, 
			:p_g1h1 => prot_g1_herb1, 
			:p_g2h2 => prot_g2_herb2,
			:fr => fr, 
			:f0 => f0, 
			:off_sd => off_sd, 
			:off_cv => off_cv,
			:sur_alt => sur_crop_alt, 
			:inv_frac => inv_frac, 
			:germ_prob => germ_prob, 
			:seed_sur => seed_sur, 
			:fec_max => fec_max, 
			:fec_dd => fec_dd, 
			:sur_spot => sur_spot, 
			:dis_rate => dis_rate, 
			:Y0 => Y0, 
			:Y_slope => Y_slope,  
			:Y_alt => Y_ALT, 
			:rep_pen => rep_pen, 
			:cost_WW => cost_WW, 
			:cost_alt => cost_ALT, 
			:cost_fal => cost_FAL, 
			:cost_plow => cost_plow, 
			:cost_spot => cost_spot, 
			:cost_herb => cost_herb_one)

	return (pop_list, reward_list, par_dict)

end












