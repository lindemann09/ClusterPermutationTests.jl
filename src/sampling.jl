function initial_fit!(cpc::CPCollection,
    data::CPData)
    @unpack specs, estimate_fnc, preprocess_fnc = cpc.def

    mtx = preprocess_fnc(data.mtx, data.design, specs)
    para = [estimate_fnc(v, data.design, specs) for v in eachcol(mtx)]
    _reset_vector!(cpc.stats, para)
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

function _do_resampling(rng::AbstractRNG,
    cpt::ClusterPermutationTest;
    n_permutations::Integer,
    progressmeter::Union{Nothing,Progress})::Vector{Vector}

    @unpack cpc, data = cpt
    @unpack specs, estimate_fnc, mass_fnc, preprocess_fnc = cpc.def
    T = eltype(cpc.stats)
    cl_ranges = cluster_ranges(cpc)
    data_mtx = data.mtx
    perm_design = copy(data.design) # shuffle always copy of design

    cl_stats_distr = Vector{T}[] # distribution of cluster-level statistics
    for _ in 1:n_permutations
        if !isnothing(progressmeter)
            next!(progressmeter)
        end
        shuffle_variable!(rng, perm_design, specs.iv)
        cl_stats = T[]
        for r in cl_ranges
            mtx = preprocess_fnc(data_mtx[:, r], perm_design, specs)
            p = [estimate_fnc(v, perm_design, specs) for v in eachcol(mtx)]
            push!(cl_stats, mass_fnc(p))
        end
        push!(cl_stats_distr, cl_stats)
    end
    return cl_stats_distr
end;

# utilities

function _reset_vector!(v::Vector, new_v::Vector)
    empty!(v)
    return append!(v, new_v)
end

;