struct CPCollection{T<:Real}
    def::CPTestDefinition # definition of test

    cc::ClusterCritODef # cluster definition

    stats::Vector{T} # test statistic of initial fit at each sample
    S::Vector{Vector{T}} # fits TODO should be matrix
end;

function CPCollection(fnc::CPTestDefinition, cluster_criterium::ClusterCritODef;
    ftype::Type{<:Float64}=Float64)

    T = ftype
    return CPCollection(fnc, cluster_criterium, T[], Vector{T}[])
end;

"""
get cluster ranges from statistics at each sample
"""
function cluster_ranges(x::CPCollection)
    return cluster_ranges(x.stats, x.cc)
end;

npermutations(x::CPCollection) = length(x.S)

function fits(x::CPCollection)::Matrix
    if length(x.S) == 0
        return zeros(eltype{x.S}, 0, 0)
    else
        return transpose(reduce(hcat, x.S))
    end
end
