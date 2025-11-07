"""
module defines fit_initial_time_series!() and resample!() for all ClusterPermutationTest

A specific test has to define

1. CP<Model> <: ClusterPermutationTest; a struct with the following fields:
	* cpc::CPCollection
	* dat::CPData
2. parameter_estimates(cpt::ClusterPermutationTest;
			fit_cluster_only::Bool=false, store_fits::Bool = false)::TParameterVector
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

function fit_initial_time_series!(cpt::ClusterPermutationTest)
	# initial fit of
	# all data samples (time_series) using (not permuted) design
	# replace existing stats
	empty!(cpt.cpc.m)
	cl_ranges = [Int32(1):Int32(epoch_length(cpt.dat))] # fit all time points
	parameter_estimates(cpt, cpt.dat.design, cl_ranges; store_fits = true)
	return nothing
end

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

	n_samples = sum(length.(cluster_ranges(cpt)))
	print("number of samples to be tested: $n_samples")
	if progressmeter
		prog = Progress(n_permutations, 0.25, "resampling")
	else
		prog = nothing
	end

	if n_threads > 1
		println(", using $n_threads threads")
		npt = convert(Int64, ceil(n_permutations/n_threads)) # n permutations per thread
		results = Vector{Vector{TParameterVector}}(undef, n_threads)
		Threads.@threads for n in 1:n_threads
			results[n] = _do_resampling(rng, cpt; n_permutations = npt, progressmeter = prog)
		end
		# combine results of all threads (each thread returns a vector of vector of parameter estimates)
		for r in results
			append!(cpt.cpc.S, r)
		end
	else
		println("")
		sampling_results = _do_resampling(rng, cpt; n_permutations, progressmeter = prog)
		append!(cpt.cpc.S, sampling_results)
	end
	return nothing
end;

function _do_resampling(rng::AbstractRNG,
	cpt::ClusterPermutationTest;
	n_permutations::Integer,
	progressmeter::Union{Nothing, Progress})::Vector{TParameterVector}

	design = copy(cpt.dat.design) # shuffle always copy of design
	cl_ranges = cluster_ranges(cpt) # get the ranges of cluster to do the permutation test #TODO console feedback?
	cl_stats_distr = TParameterVector[] # distribution of cluster-level statistics per cluster

	for _ in 1:n_permutations
		shuffle_variable!(rng, design, cpt.cpc.iv) # shuffle design
		p = parameter_estimates(cpt, design, cl_ranges) # get parameter estimates for this cluster
		push!(cl_stats_distr, p)
		isnothing(progressmeter) || next!(progressmeter)
	end
	return cl_stats_distr
end;
