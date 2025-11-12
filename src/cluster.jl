###
### Cluster definitions and criteriums
###
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

function _cluster_ranges(dat::AbstractArray{Float64}, cc::ClusterCriterium)::Vector{TClusterRange}
	# find clusters in dat according to cc
	d = cc.use_absolute ? abs.(dat) : dat
	threshold = cc.threshold
	min_size = cc.min_size
	# find cluster
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

cluster_ranges(::Any, cc::ClusterDefinition)::Vector{TClusterRange} = cc.ranges

function _cluster_mass(mass_fnc::Function, dat::AbstractArray{Float64}, cl_ranges::Vector{TClusterRange})
	# compute cluster mass for all clusters detected in dat
	return [mass_fnc(dat[cl]) for cl in cl_ranges]
end;

function _cluster_pvalues(
	perm_stats::Matrix{Float64},
	cluster_mass::Vector{Float64},
	inhibit_warning::Bool)::Vector{Float64}
	# Monte Carlo permutation p value

	n = size(perm_stats, 1)
	if n < 5_000
		if !inhibit_warning
			@warn "Small number of permutations. Estimate of p is not precise! " *
				  "Call resample!."
		end
		if n < 1_000
			return []
		end
	end

	p = []
	for (i, stats) in enumerate(cluster_mass)
		n_l = sum(perm_stats[:, i] .> abs(stats))
		append!(p, 2 * (n_l / n))
	end
	return p
end;


function _cluster_table(smpl_stats::AbstractArray{Float64},
	cl_ranges::Vector{TClusterRange},
	pvals::Vector{Float64})::Table

	if length(pvals) > 0
		p = pvals
		sign = [p <= 0.05 ? "*" : "" for p in pvals]
	else
		p = repeat(["?"], length(cl_ranges))
		sign = repeat([""], length(cl_ranges))
	end

	return Table(; id = 1:length(cl_ranges),
		from = [c.start for c in cl_ranges],
		to = [c.stop for c in cl_ranges],
		size = [c.stop - c.start + 1 for c in cl_ranges],
		min = [minimum(smpl_stats[c]) for c in cl_ranges],
		max = [maximum(smpl_stats[c]) for c in cl_ranges],
		p = p,
		sign = sign)
end


## helper function
function _join_ranges(vec::Vector{UnitRange{T}}) where T <: Integer
	rtn::Vector{UnitRange{T}} = []
	for x in sort(vec)
		if isempty(rtn) || (rtn[end].stop+1 < x.start)
			push!(rtn, x)
		else
			rtn[end] = rtn[end].start:maximum((x.stop, rtn[end].stop))
		end
	end
	return rtn
end
_join_ranges(vec::Vector{Vector{UnitRange{T}}}) where T  =
		 _join_ranges(vcat(vec...))
