"""
module defines fit_initial_time_series!() and resample!() for all ClusterPermutationTest

A specific test has to define

TODO
1. CP<Model> <: ClusterPermutationTest; a struct with the following fields:
	* cpc::CPCollection
	* dat::CPData
2. parameter_estimates(cpt::ClusterPermutationTest;
			fit_cluster_only::Bool=false, store_fits::Bool = false)::Vector{TParameterVector}

	returns vector (time) of parameter vector

	function to estimates for the entire time series for a given permutation
	data might contain different data as in cpt struct (entire time series & not permuted design),
	that is, the design might be permuted and/or epochs might be the data of merely a particular cluster.
	list of test_statistics has to be returned as TParameterVector
	if store_fits is true, the function has to store the fitted models in cpt.cpc.m
3. time_series_stats(cpt::ClusterPermutationTest)::TParameterVector
	function to extract the test statistics from the initial fit stored in cpt.cpc.m
4. StatsAPI.fit(::Type{}, ...)
	the function has to create an instance of CP<Model>, call fit_initial_time_series!(..) on it
	to detect clusters to be tested and return the instance
5. test_info(x::ClusterPermutationTest)
	function returning a string with information about the test
"""

"""initial fit of  all data samples (time_series) using (not permuted) design"""
function fit_initial_time_series!(cpt::ClusterPermutationTest)

	empty!(cpt.cpc.S)

	# replace existing fits (m) and coefs
	empty!(cpt.cpc.m)
	c = parameter_estimates(cpt, cpt.dat.design;
						fit_cluster_only = false, store_fits = true)
	cpt.cpc.coefs = stack(c, dims=1) # time X effects

	# write new time points
	all_cl_ranges = [cluster_ranges(cpt, i) for i in 1:ncoefs(cpt)]
	cpt.cpc.tp = collect(Iterators.flatten(_join_ranges(all_cl_ranges)))
	return nothing
end

#= 	FIXME MAYBE NOT MUTABLE AFTER ALL?
	# write existing fits into model fit vector cpt.cpc.m
	empty!(cpt.cpc.m) # should be empty already
	c = parameter_estimates(cpt, cpt.dat.design;
		fit_cluster_only = false, store_fits = true)

	# store coefs as matrix
	coefs = stack(c, dims = 1) # time X effects

	# write new time points
	all_cl_ranges = [cluster_ranges(cpt, i) for i in 1:ncoefs(cpt)]
	tp = collect(Iterators.flatten(_join_ranges(all_cl_ranges)))

	return CPCollection(cpt.cpc.shuffle_ivs, cpt.cpc.mass_fnc, cpt.cpc.cc,
		cpt.cpc.m, coefs, tp, TParameterMatrix[])
=#

###
### Resampling
###
resample!(cpt::ClusterPermutationTest, n_permutations::Integer; kwargs...) =
	resample!(Random.GLOBAL_RNG, cpt, n_permutations; kwargs...)

function resample!(rng::AbstractRNG,
	cpt::ClusterPermutationTest,
	n_permutations::Integer;
	progressmeter::Bool = true,
	use_threads::Union{Integer, Bool} = true)

	if use_threads === true
		n_threads = n_threads_default(cpt)
	elseif use_threads > 1 # use_threads is integer
		n_threads = min(use_threads, Threads.nthreads())
	else
		n_threads = 1
	end

	n_samples = length(cpt.cpc.tp)
	print("number of samples to be tested: $n_samples")
	if progressmeter
		prog = Progress(n_permutations, 0.25, "resampling")
	else
		prog = nothing
	end

	if n_threads > 1
		println(", using $n_threads threads")
		npt = convert(Int64, ceil(n_permutations/n_threads)) # n permutations per thread
		results = Vector{Vector{TVecTimeXParameter}}(undef, n_threads) # permutations per threads
		Threads.@threads for n in 1:n_threads
			results[n] = _do_resampling(rng, cpt; n_permutations = npt, progressmeter = prog)
		end
	else
		println("")
		results = [_do_resampling(rng, cpt; n_permutations, progressmeter = prog)]
	end

	# make matrices of each effect and store in cpt.cpc.S
	n_effects = ncoefs(cpt)
	# for each effect one vector (permutation) of vector (time) of floats)
	effect_array = [TParameterVector[] for _ in 1:n_effects]
	for thread_result in results
		for sample in thread_result
			for effect_id in 1:n_effects
				eff = getindex.(sample, effect_id) # one coefficient across time
				push!(effect_array[effect_id], eff)
			end
		end
	end

	for (cnt, effects) in enumerate(effect_array)
		mtx = stack(effects, dims=2) # make matrix (time X permutation))
		if length(cpt.cpc.S) < cnt
			push!(cpt.cpc.S, mtx)
		else
			cpt.cpc.S[cnt] = hcat(cpt.cpc.S[cnt], mtx)
		end
	end
	return nothing
end;


"""returns vector (sample) of vector (time) of vector (effect)"""
function _do_resampling(rng::AbstractRNG,
	cpt::ClusterPermutationTest;
	n_permutations::Integer,
	progressmeter::Union{Nothing, Progress})::Vector{TVecTimeXParameter}

	design = copy(cpt.dat.design) # shuffle always copy of design
	# prepare vector of effect matrix
	permutations = TVecTimeXParameter[]
	for _ in 1:n_permutations
		shuffle_variable!(rng, design, cpt.cpc.shuffle_ivs) # shuffle design
		params = parameter_estimates(cpt, design;
			fit_cluster_only = true, store_fits = false) # get parameter estimates for this cluster (time x effect)
		push!(permutations, params)
		isnothing(progressmeter) || next!(progressmeter)
	end
	return permutations
end
