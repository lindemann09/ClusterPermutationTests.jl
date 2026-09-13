const ClusterStatistics = [:clusterwise, :maxmass]

struct CPConfig
	cc::TClusterCritODef # cluster definition
	mass_fnc::Function # cluster mass function
	cluster_statistic::Symbol # clusterwise, max
    logger::Union{AbstractLogger, Nothing}
end;


"""
    CPConfig(; cluster_threshold=0, cluster_min_size=10, cluster_threshold_two_sided=true,
        predefined_cluster_definition=nothing, mass_fnc=sum, cluster_statistic=:maxmass)

    CPConfig(cc::Union{ClusterCriterium, ClusterDefinition};
            kwargs....)

Creates and returns a `CPConfig` object for cluster permutation tests.

Cluster criteria
- `cluster_threshold`: The threshold for defining clusters (default: 0).
- `cluster_min_size`: The minimum size of clusters (default: 10).
- `cluster_threshold_two_sided`: Whether to use a two-sided threshold (default: true).
- `predefined_cluster_definition`: A predefined cluster definition (default: nothing).

Either `cluster_threshold` or `predefined_cluster_definition` must be specified, but not both.
Alternatively, you can directly provide a `ClusterCriterium` or `ClusterDefinition` object.

Cluster-based statistics
- `mass_fnc`: The function to compute cluster mass (default: sum).
- `cluster_statistic`: The statistic to use for clusters (default: :maxmass).

Further options
- `logger`: An optionl logger for logging information (default: nothing).

"""
function CPConfig(cc::TClusterCritODef;
        mass_fnc::Function = sum,
        cluster_statistic::SymbolOString = :maxmass,
        logger::Union{AbstractLogger, Nothing} = nothing)

    cluster_statistic = Symbol(cluster_statistic)
    if !in(cluster_statistic, ClusterStatistics)
		throw(ArgumentError("Cluster statistic $(cluster_statistic) not supported."))
	end
    return CPConfig(cc,  mass_fnc, cluster_statistic, logger)
end

function CPConfig(;
        cluster_threshold::Real=0,
        cluster_min_size::Integer = 10,
        cluster_threshold_two_sided::Bool = true,
        predefined_cluster_definition::Union{Nothing, TClusterRange, Vector{TClusterRange}} = nothing,
        kwargs...)

    if (predefined_cluster_definition isa Nothing && cluster_threshold == 0)
        throw(ArgumentError("Please specify a cluster threshold or a predefined cluster definition."))
    elseif (!(predefined_cluster_definition isa Nothing) && cluster_threshold != 0)
        throw(ArgumentError("Please specify either a cluster threshold or a predefined cluster definition, not both."))
    elseif cluster_threshold != 0
        cc = ClusterCriterium(threshold=cluster_threshold,
                                min_size=cluster_min_size,
                                use_absolute=cluster_threshold_two_sided)
    else
        cc = ClusterDefinition(predefined_cluster_definition)
    end
    return CPConfig(cc; kwargs...)
end

