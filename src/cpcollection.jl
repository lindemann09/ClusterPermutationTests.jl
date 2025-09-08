
struct ClusterPermutationCollection{T<:Real}
    specs::NamedTuple # all parameter passed together with data to estimate_fnc

    stats::Vector{T} # test statistic at each sample
    cc::ClusterDef

    S::Vector{Vector{T}} # fits TODO should be matrix
end;

function ClusterPermutationCollection(;
    cluster_criteria::ClusterDef,
    ftype::Type{<:Float64}=Float64,
    kwargs...)
    T = ftype
    specs = (; kwargs...)
    return ClusterPermutationCollection(specs, T[], cluster_criteria, Vector{T}[])
end;

function cluster_ranges(x::ClusterPermutationCollection)
    return cluster_ranges(x.stats, x.cc)
end;

function initial_fit!(cpc::ClusterPermutationCollection;
    cpt_def::ClusterPermutationTestDefinition,
    data::CPData)
    @unpack estimate_fnc, preprocess_fnc = cpt_def
    @unpack specs = cpc

    mtx = preprocess_fnc(data.mtx, data.design, specs)
    para = [estimate_fnc(mtx[:, s], data.design, specs) for s in 1:epoch_length(data)]
    reset_vector!(cpc.stats, para)
    return nothing
end

function _resample!(rng::AbstractRNG,
    def::ClusterPermutationTestDefinition,
    cpc::ClusterPermutationCollection,
    data::CPData;
    n_permutations::Integer,
    progressmeter::Bool,
    use_threads::Bool)
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
            container[n] = _do_resampling(rng, def, cpc, data;
                n_permutations=npt, progressmeter=prog)
            prog = nothing # only first thread should have a progressbar
        end
        rstats = Vector{eltype(cpc.stats)}[] # vector of vector
        for para in container
            append!(rstats, para)
        end
    else
        rstats = _do_resampling(rng, def, cpc, data; n_permutations, progressmeter=prog)
    end
    append!(cpc.S, rstats)
    return nothing
end;

function _do_resampling(rng::AbstractRNG,
    def::ClusterPermutationTestDefinition,
    cpc::ClusterPermutationCollection,
    data::CPData;
    n_permutations::Integer,
    progressmeter::Union{Nothing,Progress})::Vector{Vector}
    T = eltype(cpc.stats)
    @unpack estimate_fnc, mass_fnc, preprocess_fnc = def
    @unpack specs = cpc
    cl_ranges = cluster_ranges(cpc)
    data_mtx = data.mtx
    perm_design = copy(data.design)

    cl_stats_distr = Vector{T}[] # distribution of cluster-level statistics
    for _ in 1:n_permutations
        if !isnothing(progressmeter)
            next!(progressmeter)
        end
        shuffle_variable!(rng, perm_design, specs.iv)
        cl_stats = T[]
        for r in cl_ranges
            mtx = preprocess_fnc(data_mtx[:, r], perm_design, specs)
            p = [estimate_fnc(Vector(v), perm_design, specs) for v in eachcol(mtx)]
            push!(cl_stats, mass_fnc(p))
        end
        push!(cl_stats_distr, cl_stats)
    end
    return cl_stats_distr
end;

npermutations(x::ClusterPermutationCollection) = length(x.S)

function fits(x::ClusterPermutationCollection)::Matrix
    if length(x.S) == 0
        return zeros(eltype{x.S}, 0, 0)
    else
        return transpose(reduce(hcat, x.S))
    end
end

# utilities

function reset_vector!(v::Vector, new_v::Vector)
    empty!(v)
    return append!(v, new_v)
end

;