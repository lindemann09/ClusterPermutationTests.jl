struct ClusterCriteria
	threshold::Real
	min_size::Int
	use_absolute::Bool
end

struct ClusterDefinition
	ranges::Vector{UnitRange}
end

ClusterDef = Union{ClusterCriteria, ClusterDefinition}

function ClusterCriteria(; threshold::Real,
	min_size::Int = 10,
	use_absolute::Bool = true)
	return ClusterCriteria(threshold, min_size, use_absolute)
end

function ClusterDefinition(single_range::UnitRange)
	return ClusterDefinition([single_range])
end


function cluster_ranges(dat::Vector{<:Real}, cc::ClusterCriteria)::Vector{UnitRange}
	d = cc.use_absolute ? abs.(dat) : dat
	@unpack threshold, min_size = cc
	## find cluster
	cluster_size = 0
	ranges = UnitRange[]
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

function cluster_ranges(::Vector{<:Real}, cc::ClusterDefinition)::Vector{UnitRange}
	return cc.ranges
end



function cluster_statistics(mass_fnc::Function,
	stats::Vector{<:Real},
	cc::ClusterDef)
	ranges = cluster_ranges(stats, cc)
	return [mass_fnc(stats[cl]) for cl in ranges]
end;
