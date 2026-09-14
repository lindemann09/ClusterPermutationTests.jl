const ClusterStatistics = ["clusterwise", "maxmass", "maxmassX"]

struct CPConfig
	cc::ClusterCriterium # cluster criterium
	mass_fnc::Function # cluster mass function
	cluster_statistic::String # clusterwise, maxmass
    mxms::Int # minimum cluster size for maxmassX
    predefined_cluster::Union{Nothing, ClusterDefinition} # predefined cluster definition
    null_logger::Bool # whether to use a null logger
end;

"""
    CPConfig(; cluster_threshold:Real,
            cluster_min_size:Integer = 10,
            cluster_threshold_two_sided:Bool = true,
            predefined_cluster::Union{Nothing, UnitRange, Vector{UnitRange}, ClusterDefinition} = nothing,
            mass_fnc::Function = sum,
            cluster_statistic::String = "maxmass")

    CPConfig(cc::Union{ClusterCriterium, ClusterDefinition};
            kwargs....)

Creates and returns a `CPConfig` object for cluster permutation tests.

Cluster criteria
- `cluster_threshold`: Real The threshold for defining clusters (Real, default: 0).
- `cluster_min_size`: The minimum size of clusters (Integer, default: 10).
- `cluster_threshold_two_sided`: Whether to use a two-sided threshold (Bool,default: true).
- `predefined_cluster`: A optional predefined cluster definition
            (UnitRange, Vector{UnitRange}, ClusterDefinition).

Either `cluster_threshold` or `predefined_cluster` must be specified, but not both.
Alternatively, you can directly provide a `ClusterCriterium` or `ClusterDefinition` object.

Cluster-based statistics
- `mass_fnc`: The function to compute cluster mass (default: sum).
- `cluster_statistic`: The statistic to use for clusters in the permutation test(String, default: "maxmass").
  - "clusterwise": Cluster mass of each cluster interval after the permutation.
  - "maxmass": Maximum of cluster mass of all clusters in this permutation.
  - "maxmassX": as above, but consider on cluster with a minimum size of `X` (where X must be a positive integer).

Further options
- `null_logger`: A boolean indicating whether to use a null logger (default: true).

"""
function CPConfig(
        cc::ClusterCriterium;
        mass_fnc::Function = sum,
        cluster_statistic::SymbolOString = "maxmass",
        predefined_cluster::Union{Nothing, TClusterRange, Vector{TClusterRange}, ClusterDefinition} = nothing,
        null_logger::Bool = true)

    cluster_statistic = String(cluster_statistic)
    if !(predefined_cluster == nothing || predefined_cluster isa ClusterDefinition)
        predefined_cluster = ClusterDefinition(predefined_cluster)
    end

    # find minimum cluster size for maxmassX, or zero if not a maxmass label
    mxms = 0
    if is_maxmass(cluster_statistic)
        min_size_str = replace(cluster_statistic, "maxmass" => "")
        if isempty(min_size_str)
            mxms = 2
        else
            try
                mxms = parse(Int, min_size_str)
            catch e
                throw(ArgumentError("Invalid maxmassX label: $label. X must be a positive integer and not $min_size_str."))
            end
            if mxms < 2
                throw(ArgumentError("min_size must be an integer > 1"))
            end
        end
    elseif !in(cluster_statistic, ClusterStatistics)
        # not a valid other method
		throw(ArgumentError("Cluster statistic $(cluster_statistic) not supported. Please choose one of $(ClusterStatistics)."))
	end
    return CPConfig(cc,  mass_fnc, cluster_statistic, mxms, predefined_cluster, null_logger)
end

function CPConfig(;
        cluster_threshold::Real,
        cluster_min_size::Integer = 2,
        cluster_threshold_two_sided::Bool = true,
        kwargs...)
    cc = ClusterCriterium(threshold=cluster_threshold,
                        min_size=cluster_min_size,
                        use_absolute=cluster_threshold_two_sided)
    return CPConfig(cc; kwargs...)
end

function cluster_type(cp_config::CPConfig)::TClusterCritODef
    if cp_config.predefined_cluster !== nothing
		return cp_config.predefined_cluster
	else
		return cp_config.cc
	end
end

is_maxmass(s::String)::Bool = startswith(s, "maxmass")
is_maxmass(cp_config::CPConfig)::Bool = is_maxmass(cp_config.cluster_statistic)
