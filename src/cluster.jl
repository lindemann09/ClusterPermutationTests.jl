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

function _clusterranges(dat::AbstractArray{Float64}, cc::ClusterCriterium)::Vector{TClusterRange}
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
				# valid cluster until the end (as if cluster ends c+1)
				push!(ranges, (c+1-cluster_size):c)
			end
		end
	end
	return ranges
end;

_clusterranges(::Any, cc::ClusterDefinition)::Vector{TClusterRange} = cc.ranges

function _cluster_mass_stats(mass_fnc::Function, dat::AbstractArray{Float64}, cl_ranges::Vector{TClusterRange})
	# compute cluster mass for all clusters detected in dat
	return [mass_fnc(dat[cl]) for cl in cl_ranges]
end;

function _cluster_pvalues(
	cl_nhd::Matrix{Float64},
	cl_mass_stats::Vector{Float64},
	inhibit_warning::Bool;
	one_tail::Bool = false)::Vector{Float64}
	# Monte Carlo permutation p value

	n = size(cl_nhd, 1)
	if n < 5_000
		if !inhibit_warning
			@warn "Small number of permutations. Estimate of p is not precise! " *
				  "Call resample!."
		end
		if n < 1_000
			return []
		end
	end

	rtn = []
	for (nhd, cms) in zip(eachcol(cl_nhd), cl_mass_stats)
		p = 1 - quantilerank(abs.(nhd), abs(cms); method = :exc)
		if one_tail
			p = p / 2
		end
		append!(rtn, p)
	end
	return rtn
end;



function _cluster_table(coef_id::Integer,
	coefname::String,
	cl_ranges::Vector{TClusterRange},
	cms::Vector{Float64},
	pvals::Vector{Float64};
	add_effect_names::Bool = false,
	sign_level::Real = 0.05)::CoefTable

	if length(pvals) > 0
		p = pvals
		sign = [p <= sign_level ? "*" : "" for p in pvals]
	else
		p = repeat([NaN], length(cl_ranges))
		sign = repeat([""], length(cl_ranges))
	end
	cid = collect(1:length(cl_ranges))
	#names = repeat([coefname], length(cl_ranges))
	from = [c.start for c in cl_ranges]
	to = [c.stop for c in cl_ranges]
	size = [c.stop - c.start + 1 for c in cl_ranges]
	colnms = ["cluster", "from", "to", "size", "mass stats", "Pr(>|z|)", "sign"]
	cols = [cid, from, to, size, cms, p, sign]
	rownms = [string(coef_id) * "."*string(i) for i in cid]
	if add_effect_names && !isempty(rownms)
		rownms[1] *= " - " * coefname
	end
	pvalcol = 6
	teststatcol = 5

	return CoefTable(cols, colnms, rownms, pvalcol, teststatcol)
end

function _joined_ranges(vec::Vector{TClusterRange})
	rtn::Vector{TClusterRange} = []
	for x in sort(vec)
		if isempty(rtn) || (rtn[end].stop+1 < x.start)
			push!(rtn, x)
		else
			rtn[end] = rtn[end].start:maximum((x.stop, rtn[end].stop))
		end
	end
	return collect(Iterators.flatten(rtn))
end
_joined_ranges(vec::Vector{Vector{TClusterRange}}) = _joined_ranges(vcat(vec...))
