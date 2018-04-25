# some helper function to parse and unpack the solver results

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

# takes an array of rewards, each element and array of rewards at each
# generation. Return the indicies of n best seqences in each generation
function best_n_seq(n::Int64, reward_arr::Array{Array, 1})

	N_gen = length(reward_arr)
	best_ind = Array{Int64, 2}(N_gen, n)

	for g in 1:N_gen

		si = sortperm(reward_arr[g], rev = true)
		best_ind[g, :] = si[1:n]

	end

	return best_ind

end

# pull out the rewards
function reward_gen(reward_arr::Array{Array, 1}, ind::Array{Int64, 2})

	N_gen = length(reward_arr)
	N_series = size(ind)[2]

	out = zeros(N_gen, N_series)
	# run through each series by generation, pulling out the rewards
	for g in 1:N_gen

		out[g, :] = reward_arr[g][ind[g, :]]

	end

	return out

end

# turn a seqence of action mumbers to seqences of herb, crop, plow, and spot
function actseq_subact(actseq::Array{Int64, 2}, A::Tuple)

	# holding arrays
	herb_seq = Array{Int64, 2}(size(actseq))
	crop_seq = Array{Int64, 2}(size(actseq))
	plow_seq = Array{Int64, 2}(size(actseq))
	spot_seq = Array{Int64, 2}(size(actseq))

	for s in 1:size(actseq)[1]

		acts = act_seq_2_sub_act(A, actseq[s, :])

		herb_seq[s, :] = acts[:, ACT_HERB]
		crop_seq[s, :] = acts[:, ACT_CROP]
		plow_seq[s, :] = acts[:, ACT_PLOW]
		spot_seq[s, :] = acts[:, ACT_SPOT]

	end

	return (herb_seq, crop_seq, plow_seq, spot_seq)

end

# take the whole solution and output the best performing seqeunce in each
# generation, split into herb, crop, plow and spot
function get_best_seq(sol_ob::Tuple, A::Tuple)

	best_ind = best_n_seq(1, sol_ob[2])

	N_gen = length(sol_ob[1])
	N_t = size(sol_ob[1][1])[2]

	herb = Array{Int64, 2}(N_gen, N_t)
	crop = Array{Int64, 2}(N_gen, N_t)
	plow = Array{Int64, 2}(N_gen, N_t)
	spot = Array{Int64, 2}(N_gen, N_t)

	for g in 1:N_gen

		acts = act_seq_2_sub_act(A, sol_ob[1][g][best_ind[g], :])

		herb[g, :] = acts[:, ACT_HERB]
		crop[g, :] = acts[:, ACT_CROP]
		plow[g, :] = acts[:, ACT_PLOW]
		spot[g, :] = acts[:, ACT_SPOT]

	end

	return (herb, crop, plow, spot)

end

# show best single sequence, formatted for compact display using ggplot2
# in R for easier faceting and labelling options (longform data.frame)
function best_seq_2_df(sol_ob::Tuple, A::Tuple)

	best_ind = best_n_seq(1, sol_ob[2])

	N_t = size(sol_ob[1][1])[2]

	acts = act_seq_2_sub_act(A, sol_ob[1][end][best_ind[end], :])

	# make into long form dataframe for plotting
	return DataFrame(sub_act = vcat(repmat(["herbicide"], N_t),
			repmat(["crop"], N_t), repmat(["plow"], N_t),
			repmat(["spot"], N_t)),
		ts = repmat(1:N_t, 4),
		best_act = vcat(acts[:, ACT_HERB],
			acts[:, ACT_CROP], acts[:, ACT_PLOW],
			acts[:, ACT_SPOT]))

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

function sub_act_2_act_seq(A::Tuple, sub_act::Array{Int64, 2})

	# put the Action space into a vectorised form to help searching
	A_mat = reshape(vcat(A...), 4, length(A))

	act_seq = Array{Int64, 1}(T)

	for t in 1:T

		herb_act = find(A_mat[ACT_HERB, :] .== sub_act[t, ACT_HERB])
		crop_act = find(A_mat[ACT_CROP, :] .== sub_act[t, ACT_CROP])
		plow_act = find(A_mat[ACT_PLOW, :] .== sub_act[t, ACT_PLOW])
		spot_act = find(A_mat[ACT_SPOT, :] .== sub_act[t, ACT_SPOT])

		act_seq[t] = intersect(herb_act, crop_act, plow_act, spot_act)[1]

	end

	return act_seq

end
