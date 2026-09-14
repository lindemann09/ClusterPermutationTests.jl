"""
module defines fit_initial_time_series!() and resample!() for all ClusterPermutationTest

A specific test has to define

TODO
1. CP<Model> <: ClusterPermutationTest; a struct with the following mandatory fields:
	* config::CPConfig
	* cpc::CPCollection
	* dat::CPData

2. parameter_estimates(cpt::ClusterPermutationTest,
			design::AbstractStudyDesign,
			time_points::Vector{<:Integer};
			store_model_fits::Bool = false)::T2DParamVector

	returns vector time x parameter

	function to estimates for the entire time series for a given permutation
	data might contain different data as in cpt struct (entire time series & not permuted design),
	that is, the design might be permuted and/or epochs might be the data of merely a particular cluster.
	list of test_statistics has to be returned as TParameterVector
	if store_model_fits is true, the function has to store the fitted models in cpt.cpc.M

3. time_series_stats(cpt::ClusterPermutationTest)::TParameterVector
	function to extract the test statistics from the initial fit stored in cpt.cpc.M

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
	empty!(cpt.cpc.M)

	old_logger = cpt.config.null_logger ? global_logger(NullLogger()) : nothing # change logger

	atp = collect(1:epoch_length(cpt.dat)) # all time points
	c = parameter_estimates(cpt, cpt.dat.design, atp; store_model_fits = true)

	isnothing(old_logger) || global_logger(old_logger)

	cpt.cpc.coefs = stack(c, dims = 1) # time X effects
	return nothing
end


###
### Resampling
###
"""
    resample!(cpt::ClusterPermutationTest, n_permutations::Integer; kwargs...)
    resample!(rng::AbstractRNG, cpt::ClusterPermutationTest, n_permutations::Integer;
              progressmeter=nothing, use_threads=true, logger=NullLogger())

Run `n_permutations` permutations and accumulate the null-hypothesis distribution in `cpt`.

Repeated calls append to the existing permutation samples.

# Keyword arguments
- `progressmeter`: show a progress bar; defaults to `true` when running in a terminal.
- `use_threads`: `true` uses all available threads, an integer limits the count, `false`/`1`
  uses a single thread.
"""
resample!(cpt::ClusterPermutationTest, n_permutations::Integer; kwargs...) =
	resample!(Random.GLOBAL_RNG, cpt, n_permutations; kwargs...)

function resample!(rng::AbstractRNG,
	cpt::ClusterPermutationTest,
	n_permutations::Integer;
	progressmeter::Union{Bool, Nothing} = nothing,
	use_threads::Union{Integer, Bool} = true)
	# progressmeter used per default, if terminal is used (not jupyter or quarto)

	if use_threads === true
		n_threads = n_threads_default(cpt)
	elseif use_threads > 1 # use_threads is integer
		n_threads = min(use_threads, Threads.nthreads())
	else
		n_threads = 1
	end

	if cpt.config.cluster_statistic == "clusterwise"
		all_cluster = _cluster_ranges(cpt.cpc.coefs, cluster_type(cpt.config))
		n_samples = length(_joined_ranges(all_cluster))
	else
		n_samples = epoch_length(cpt.dat)
	end
	print("To-be tested samples ($(cpt.config.cluster_statistic)): $n_samples")


	if progressmeter === nothing
		progressmeter = isa(stderr, Base.TTY) # if terminal is used and not jupyter or Quarto
	end
	if progressmeter
		prog = Progress(n_permutations; dt = 0.25, desc = "Resampling")
	else
		prog = nothing
	end

	old_logger = cpt.config.null_logger ? global_logger(NullLogger()) : nothing # change logger

	# result:  is a Vector thread x permutation x effect x cluster
	if n_threads > 1
		println(", using $n_threads threads")
		npt = convert(Int64, ceil(n_permutations/n_threads)) # n permutations per thread
		if cpt.config.cluster_statistic == "clusterwise"
			results = Vector{Vector{T2DParamVector}}(undef, n_threads) # permutations per threads
			Threads.@threads for n in 1:n_threads
				results[n] = _resampling_cluster_wise(rng, cpt, npt, prog)
			end
		else
			results = Vector{T2DParamVector}(undef, n_threads) # permutations per threads
			Threads.@threads for n in 1:n_threads
				results[n] = _resampling_max_cluster_stats(rng, cpt, npt, prog)
			end
		end

		Threads.@threads for n in 1:n_threads
			if cpt.config.cluster_statistic == "clusterwise"
				results[n] = _resampling_cluster_wise(rng, cpt, npt, prog)
			else
				results[n] = _resampling_max_cluster_stats(rng, cpt, npt, prog)
			end
		end
	else
		println("")
		if cpt.config.cluster_statistic == "clusterwise"
			results = [_resampling_cluster_wise(rng, cpt, n_permutations, prog)]
		else
			results = [_resampling_max_cluster_stats(rng, cpt, n_permutations, prog)]
		end
	end
	isnothing(old_logger) || global_logger(old_logger)
	isnothing(prog) || finish!(prog)

	n_effects = ncoefs(cpt)
	effects_cl_masses = [T2DParamVector() for _ in 1:n_effects]

	## combine all threads results
	# make 3d-vector of cluster masses: effect x all combined permutations x clusters (n cluster=1 for maxmass)
	for thread_result in results
		for permutation in thread_result
			for (eid, cms_eff) in enumerate(permutation)
				if cpt.config.cluster_statistic == "clusterwise"
					push!(effects_cl_masses[eid], cms_eff)
				else
					push!(effects_cl_masses[eid], [cms_eff])
				end
			end
		end
	end

	# make one matrix per effect and store in cpt.cpc.S
	# append cpt.cpc.S:  vector (effect) of matrix sample X cluster (with cluster=1 for max)
	append_samples = length(cpt.cpc.S) > 0
	for (eid, effects) in enumerate(effects_cl_masses)
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
@inline function _resampling_max_cluster_stats(rng::AbstractRNG,
	cpt::ClusterPermutationTest,
	n_permutations::Integer,
	progressmeter::Union{Nothing, Progress})::T2DParamVector # vector (permutations) x effect

	design = copy(cpt.dat.design) # shuffle always copy of design

	time_points = collect(1:epoch_length(cpt.dat)) # all time points
	n_effects = ncoefs(cpt)

	cc = ClusterCriterium(
		threshold = cpt.config.cc.threshold,
		min_size = cpt.config.mxms, # different min_size for permutations
		use_absolute = cpt.config.cc.use_absolute)

	# prepare vector (permutation) x effect
	permutations = TParameterVector[]
	for _ in 1:n_permutations
		shuffle_variable!(rng, design, cpt.cpc.shuffle_ivs) # shuffle design
		# get parameter estimates for the time points (time x effect)
		params = parameter_estimates(cpt, design, time_points; store_model_fits = false)

		max_cluster_masses = TParameterVector(undef, n_effects)
		for eid in 1:n_effects
			ts = getindex.(params, eid) # time series stats for this effect
			cms = _cluster_mass_stats(cpt.config.mass_fnc, ts, _cluster_ranges(ts, cc))
			max_cluster_masses[eid] = isempty(cms) ? 0.0 : maximum(cms)
		end
		push!(permutations, max_cluster_masses)

		isnothing(progressmeter) || next!(progressmeter)
	end
	return permutations
end

"""returns vector (sample) of vector (time) of vector (effect)"""
@inline function _resampling_cluster_wise(rng::AbstractRNG,
	cpt::ClusterPermutationTest,
	n_permutations::Integer,
	progressmeter::Union{Nothing, Progress})::Vector{T2DParamVector}
	# vector (permutations) x effect X cluster

	design = copy(cpt.dat.design) # shuffle always copy of design

	# prepare vector (cms) of effect x cluster
	permutations = T2DParamVector[]

	all_cluster = _cluster_ranges(cpt.cpc.coefs, cluster_type(cpt.config))
	time_points = _joined_ranges(all_cluster)

	# idx: ranges of indices for the returns parameters that correspond to the time points in the cluster
	idx = deepcopy(all_cluster)  # allocated memory
	for i in eachindex(idx)
		for j in eachindex(idx[i])
			cl = all_cluster[i][j]
			idx[i][j] = findfirst(isequal(cl.start), time_points):findfirst(isequal(cl.stop), time_points)
		end
	end

	for _ in 1:n_permutations
		shuffle_variable!(rng, design, cpt.cpc.shuffle_ivs) # shuffle design
		# get parameter estimates for the time points (time x effect)
		params = parameter_estimates(cpt, design, time_points; store_model_fits = false)

		cms_vec = T2DParamVector()
		for (eid, effect_cluster) in enumerate(idx)
			# cluster mass statistics for each cluster of this effect
			cms = [cpt.config.mass_fnc(getindex.(params[cl_idx], eid)) for cl_idx in effect_cluster]
			push!(cms_vec, cms)
		end
		push!(permutations, cms_vec)
		isnothing(progressmeter) || next!(progressmeter)
	end
	return permutations
end
