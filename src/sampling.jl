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


function fit_initial_time_series!(cpt::ClusterPermutationTest;
		logger::Union{AbstractLogger, Nothing} = nothing)
	empty!(cpt.cpc.S)

	# replace existing fits (m) and coefs
	empty!(cpt.cpc.m)
	atp = collect(Int32(1):Int32(epoch_length(cpt.dat))) # all time points

	old_logger = isnothing(logger) ? nothing : global_logger(logger) # change logger
	c = parameter_estimates(cpt, cpt.dat.design, atp; store_fits = true)
	isnothing(old_logger) || global_logger(old_logger)

	cpt.cpc.coefs = stack(c, dims = 1) # time X effects


	# write new time points
	cpt.cpc.cl = [_cluster_ranges(d, cpt.cpc.cc) for d in eachcol(cpt.cpc.coefs)]
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
	use_threads::Union{Integer, Bool} = true,
	logger::Union{AbstractLogger, Nothing} = NullLogger())

	if use_threads === true
		n_threads = n_threads_default(cpt)
	elseif use_threads > 1 # use_threads is integer
		n_threads = min(use_threads, Threads.nthreads())
	else
		n_threads = 1
	end

	n_samples = length(_joined_ranges(cpt.cpc.cl))
	print("number of samples to be tested: $n_samples")

	if progressmeter
		prog = Progress(n_permutations, 0.25, "resampling")
	else
		prog = nothing
	end

	old_logger = isnothing(logger) ? nothing : global_logger(logger) # change logger

	if n_threads > 1
		println(", using $n_threads threads")
		npt = convert(Int64, ceil(n_permutations/n_threads)) # n permutations per thread
		results = Vector{Vector{T2DParamVector}}(undef, n_threads) # permutations per threads
		Threads.@threads for n in 1:n_threads
			results[n] = _do_resampling(rng, cpt, npt, prog)
		end
	else
		println("")
		results = [_do_resampling(rng, cpt, n_permutations, prog)]
	end
	isnothing(old_logger) || global_logger(old_logger)
	isnothing(prog) || finish!(prog)

	# make matrices of each effect and store in cpt.cpc.S
	n_effects = ncoefs(cpt)

	# thread_result is a NDVector thread x samples x effect x cluster

	# effect_array: NDVector effect x all sample X cluster
	effect_array = [T2DParamVector() for _ in 1:n_effects]
	for thread_result in results
		for sample in thread_result
			for (eid, cms_eff) in enumerate(sample)
				push!(effect_array[eid], cms_eff)
			end
		end
	end

	# append cpt.cpc.S:  vector (effect) of matrix sample X cluster
	append_samples = length(cpt.cpc.S) > 0
	for (eid, effects) in enumerate(effect_array)
		mtx = stack(effects, dims = 1) # make matrix (sample X permutation)
		if append_samples
			cpt.cpc.S[eid] = vcat(cpt.cpc.S[eid], mtx) # append it parameter exist
		else
			push!(cpt.cpc.S, mtx)
		end
	end
	return nothing
end;


"""returns vector (sample) of vector (time) of vector (effect)"""
function _do_resampling(rng::AbstractRNG,
	cpt::ClusterPermutationTest,
	n_permutations::Integer,
	progressmeter::Union{Nothing, Progress})::Vector{T2DParamVector} # vector of effect X cluster

 	design = copy(cpt.dat.design) # shuffle always copy of design

	# prepare vector (cms) of effect x cluster
	permutations = T2DParamVector[]

	time_points = _joined_ranges(cpt.cpc.cl)
	# idx: ranges of indices for the returns parameters that correspond to the time points in the cluster
	idx = deepcopy(cpt.cpc.cl)  # allocated memory
	for i in eachindex(idx)
		for j in eachindex(idx[i])
			cl = cpt.cpc.cl[i][j]
			idx[i][j] = findfirst(isequal(cl.start), time_points):findfirst(isequal(cl.stop), time_points)
		end
	end

	for _ in 1:n_permutations
		shuffle_variable!(rng, design, cpt.cpc.shuffle_ivs) # shuffle design
		# get parameter estimates for the time points (time x effect)
		params_time = parameter_estimates(cpt, design, time_points; store_fits = false)

		cms_vec = T2DParamVector()
		for (eid, effect_cluster) in enumerate(idx)
			# cluster mass statistics for each cluster of this effect
			cms = [cpt.cpc.mass_fnc(getindex.(params_time[cl_idx], eid)) for cl_idx in effect_cluster]
			push!(cms_vec, cms)
		end
		push!(permutations, cms_vec)
		isnothing(progressmeter) || next!(progressmeter)
	end
	return permutations
end
