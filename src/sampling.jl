function initial_fit!(cpt::ClusterPermutationTest)
    # initial fit of
    # all data samples (time_series) using (not permuted) design
    para = _parameter_estimate(cpt, cpt.dat.mtx, cpt.dat.design)
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

@inline function _parameter_estimate(cpt::ClusterPermutationTest, mtx::Matrix{<:Real},
	                        permutation::PermutationDesign)
    # Estimate parameters for a single time series sample
    mtx, design_tbl = prepare_data(cpt, mtx, permutation)
    return [estimate(cpt, s, design_tbl) for s in eachcol(mtx)]
end

@inline function _do_resampling(rng::AbstractRNG,
    cpt::ClusterPermutationTest;
    n_permutations::Integer,
    progressmeter::Union{Nothing,Progress})::Vector{Vector}

    @unpack cpc, dat = cpt
    mass_fnc = cpc.mass_fnc
    T = eltype(cpc.stats)
    cl_ranges = cluster_ranges(cpc)
    perm_design = copy(dat.design) # shuffle always copy of design

    cl_stats_distr = Vector{T}[] # distribution of cluster-level statistics
    for _ in 1:n_permutations
        if !isnothing(progressmeter)
            next!(progressmeter)
        end
        shuffle_variable!(rng, perm_design, cpt.iv)
        cl_stats = T[]
        for r in cl_ranges
            p = _parameter_estimate(cpt, dat.mtx[:, r], perm_design) # FIXME VIEW?
            push!(cl_stats, mass_fnc(p))
        end
        push!(cl_stats_distr, cl_stats)
    end
    return cl_stats_distr
end;
