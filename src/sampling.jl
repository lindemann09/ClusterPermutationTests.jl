"""
module defines initial_fit!() and resample!() for all ClusterPermutationTest

A specific test has to define

1. CP<Model> <: ClusterPermutationTest; a struct with the following fields:
	* cpc::CPCollection
	* dat::CPData
	* m::Vector{StatisticalModel} # fitted models of initial fit

2. parameter_estimates(cpt::ClusterPermutationTest, dat::CPData; initial_fit::Bool = false)::TParameterVector
	function to estimates for the entire time series for a given permutation
	data might contain different data as is cpt (entire time series & not permuted design),
	that is, the design might be permuted and mtx might be the data of merely a particular cluster

	list of test_statistics has to be returned as TParameterVector
	if initial_fit is true,
		* the function has to store the fitted models in cpt.cpc.m
		* the function has to store the test statistics in cpt.cpc.stats

3. StatsAPI.fit(::Type{}, ...)
	the function has to create an instance of CP<Model>, call initial_fit!(..) on it
	to detect clusters to be tested and return the instance
"""

function initial_fit!(cpt::ClusterPermutationTest)
	# initial fit of
	# all data samples (time_series) using (not permuted) design
	# replace existing stats
	empty!(cpt.cpc.stats)
	empty!(cpt.cpc.m)
	parameter_estimates(cpt, cpt.dat; initial_fit = true)
	return nothing
end

###
### Resampling
###
resample!(cpt::ClusterPermutationTest; n_permutations::Integer, kwargs...) =
	resample!(Random.GLOBAL_RNG, cpt, n_permutations; kwargs...)
resample!(rng::AbstractRNG, cpt::ClusterPermutationTest; n_permutations::Integer, kwargs...) =
	resample!(rng, cpt, n_permutations; kwargs...)
resample!(cpt::ClusterPermutationTest, n_permutations::Integer; kwargs...) =
	resample!(Random.GLOBAL_RNG, cpt, n_permutations; kwargs...)

function resample!(rng::AbstractRNG,
	cpt::ClusterPermutationTest,
	n_permutations::Integer;
	progressmeter::Bool = true,
	use_threads::Bool = true)

	if progressmeter
		n_cluster = length(cluster_ranges(cpt.cpc))
		prog = Progress(n_permutations * n_cluster, 0.25, "resampling")
	else
		prog = nothing
	end
	if use_threads
		n_thr = Threads.nthreads()
		npt = convert(Int64, ceil(n_permutations/n_thr)) # n permutations per thread
		results = Vector{Vector{TParameterVector}}(undef, n_thr)
		Threads.@threads for n in 1:n_thr
			results[n] = _do_resampling(rng, cpt; n_permutations = npt, progressmeter = prog)
		end
		# combine results of all threads (each thread returns a vector of vector of parameter estimates)
		for r in results
			_append_sampling_results!(cpt.cpc.S, r)
		end
	else
		sampling_results = _do_resampling(rng, cpt; n_permutations, progressmeter = prog)
    	_append_sampling_results!(cpt.cpc.S, sampling_results)
	end
	return nothing
end;

function _do_resampling(rng::AbstractRNG,
	cpt::ClusterPermutationTest;
	n_permutations::Integer,
	progressmeter::Union{Nothing, Progress})::Vector{TParameterVector}

	mass_fnc = cpt.cpc.mass_fnc
	design = copy(cpt.dat.design) # shuffle always copy of design
	cl_ranges = cluster_ranges(cpt) # get the ranges of cluster to do the permutation test #TODO console feedback?
	cl_stats_distr = TParameterVector[] # distribution of cluster-level statistics per cluster

	for (i, r) in enumerate(cl_ranges)
		# cluster data with shuffled design
		dat = CPData(cpt.dat.mtx[:, r], design) # TODO view?
		push!(cl_stats_distr, TParameterVector())
		for _ in 1:n_permutations
			isnothing(progressmeter) || next!(progressmeter)
			shuffle_variable!(rng, dat.design, cpt.cpc.iv) # shuffle design
			p = parameter_estimates(cpt, dat) # get parameter estimates for this cluster
			push!(cl_stats_distr[i], mass_fnc(p))
		end
	end
	return cl_stats_distr
end;

function _append_sampling_results!(a::Vector{TParameterVector}, b::Vector{TParameterVector})
	if isempty(a)
		for _ in 1:length(b)
			push!(a, TParameterVector())
		end
	end
	for (x, y) in zip(a, b)
		append!(x, y)
	end
	return nothing
end;