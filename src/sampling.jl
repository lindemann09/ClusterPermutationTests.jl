"""
module defines initial_fit!() and resample!() for all ClusterPermutationTest

A specific test has to define

1. CP<Model> <: ClusterPermutationTest; a struct with the following fields:
	* cpc::CPCollection
	* dat::CPData
	* iv::Symbol

2. parameter_estimates(cpt::ClusterPermutationTest, dat::CPData)::TParameterVector
    function to estimates for the entire time series
3. parameter_estimates(cpt::ClusterPermutationTest, dat::CPData, range::TClusterRange)::TParameterVector
    estimates for a specific section in the time series (cluster)

4. StatsAPI.fit(::Type{}, ...)
    the function has to create an instance of CP<Model>, call initial_fit!(..) on it
    to detect clusters to be tested and return the instance

Notes:
* TParameterVector is a Vector{Float64}
* TClusterRange is a UnitRange{Int32}

"""

function initial_fit!(cpt::ClusterPermutationTest)
    # initial fit of
    # all data samples (time_series) using (not permuted) design
    para = parameter_estimates(cpt, cpt.dat)
    # replace existing stats
    empty!(cpt.cpc.stats)
    append!(cpt.cpc.stats, para)
    return nothing
end

resample!(cpt::ClusterPermutationTest; kwargs...) = resample!(Random.GLOBAL_RNG, cpt; kwargs...)

function resample!(rng::AbstractRNG,
    cpt::ClusterPermutationTest;
	n_permutations::Integer,
	progressmeter::Bool = true,
	use_threads::Bool = true)

    if progressmeter
        prog = Progress(n_permutations, 1, "resampling")
    else
        prog = nothing
    end
    if use_threads
        n_thr = Threads.nthreads()
        npt = convert(Int64, ceil(n_permutations/n_thr)) # n permutations per thread
        container = Vector{Any}(undef, n_thr)
        Threads.@threads for n in 1:n_thr
            container[n] = _do_resampling(rng, cpt;
                n_permutations=npt, progressmeter=prog)
            prog = nothing # only first thread should have a progressbar
        end
        rstats = Vector{eltype(cpt.cpc.stats)}[] # vector of vector
        for para in container
            append!(rstats, para)
        end
    else
        rstats = _do_resampling(rng, cpt; n_permutations, progressmeter=prog)
    end
    append!(cpt.cpc.S, rstats)
    return nothing
end;

@inline function _do_resampling(rng::AbstractRNG,
    cpt::ClusterPermutationTest;
    n_permutations::Integer,
    progressmeter::Union{Nothing,Progress})::Vector{Vector}

    mass_fnc = cpt.cpc.mass_fnc
    T = eltype(cpt.cpc.stats)
    cl_ranges = cluster_ranges(cpt) # get the ranges of cluster to do the permutation test #TODO console feedback?

    dat = CPData(cpt.dat.mtx, copy(cpt.dat.design)) # shuffle always copy of design
    cl_stats_distr = Vector{T}[] # distribution of cluster-level statistics

    for _ in 1:n_permutations
        if !isnothing(progressmeter)
            next!(progressmeter)
        end
        shuffle_variable!(rng, dat.design, cpt.iv)
        cl_stats = T[]
        for r in cl_ranges # loop over cluster
            p = parameter_estimates(cpt, dat, r) # FIXME VIEW?
            push!(cl_stats, mass_fnc(p))
        end
        push!(cl_stats_distr, cl_stats)
    end

    return cl_stats_distr
end;
