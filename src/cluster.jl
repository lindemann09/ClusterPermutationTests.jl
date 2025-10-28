const TParameterVector = Vector{Float64}
const TClusterRange = UnitRange{Int32}

struct ClusterCriterium
	threshold::Real
	min_size::Int32
	use_absolute::Bool
end

struct ClusterDefinition
	ranges::Vector{TClusterRange}
end

const TClusterCritODef = Union{ClusterCriterium, ClusterDefinition}

function ClusterCriterium(; threshold::Real,
	min_size::Int = 10,
	use_absolute::Bool = true)
	return ClusterCriterium(threshold, min_size, use_absolute)
end

function ClusterDefinition(single_range::UnitRange)
	return ClusterDefinition([single_range])
end


function cluster_ranges(dat::TParameterVector, cc::ClusterCriterium)::Vector{TClusterRange}
	# find clusters in dat according to cc
	d = cc.use_absolute ? abs.(dat) : dat
	@unpack threshold, min_size = cc
	## find cluster
	cluster_size = 0
	ranges = TClusterRange[]
	n = length(d)
	for (c, x) in enumerate(d)
		if x >= threshold
			cluster_size += 1
		else
			if cluster_size >= min_size
				# valid new cluster
				push!(ranges, (c-cluster_size):(c-1))
			end
			cluster_size = 0
		end

		if c == n
			# last element processed
			if cluster_size >= min_size
				# valid cluster at the end
				push!(ranges, (c-cluster_size):c)
			end
		end
	end
	return ranges
end;

cluster_ranges(::Any, cc::ClusterDefinition)::Vector{TClusterRange} =
		cc.ranges

function cluster_stats(mass_fnc::Function,
	stats::TParameterVector,
	cc::TClusterCritODef)
	ranges = cluster_ranges(stats, cc)
	return [mass_fnc(stats[cl]) for cl in ranges]
end;
