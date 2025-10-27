###
### CPCollection
###

struct CPCollection{M}
	iv::Symbol # name of the to be shuffled independent variable
	mass_fnc::Function # cluster mass function
	cc::TClusterCritODef # cluster definition

	stats::TParameterVector # test statistic of initial fit at each sample
	m::Vector{M} # fitted models of initial fit
	S::Vector{TParameterVector} # samples
end;

CPCollection{M}(iv::SymbolOString, mass_fnc::Function, cluster_criterium::TClusterCritODef) where {M} =
	CPCollection(Symbol(iv), mass_fnc, cluster_criterium, TParameterVector(),
			M[], TParameterVector[])


"""
get cluster ranges from statistics at each sample
"""
cluster_ranges(x::CPCollection) = cluster_ranges(x.stats, x.cc)

npermutations(x::CPCollection) = length(x.S) == 0 ? 0 : length(x.S[1])

function fits(x::CPCollection)::Matrix
	if length(x.S) == 0
		return zeros(eltype(TParameterVector), 0, 0)
	else
		return reduce(hcat, x.S)
	end
end

###
### ClusterPermutationTest
###

abstract type ClusterPermutationTest end
#requires
#	cpc::CPCollection
#	dat::CPData
#	iv::Symbol

nepochs(x::ClusterPermutationTest) = nepochs(x.dat)
epoch_length(x::ClusterPermutationTest) = epoch_length(x.dat)
design_table(x::ClusterPermutationTest) = design_table(x.dat)

npermutations(x::ClusterPermutationTest) = npermutations(x.cpc)
StudyDesigns.unit_observation(x::ClusterPermutationTest) = unit_observation(x.dat.design.uo)
cluster_ranges(x::ClusterPermutationTest) = cluster_ranges(x.cpc)
cluster_criterium(x::ClusterPermutationTest) = x.cpc.cc
fits(x::ClusterPermutationTest) = fits(x.cpc)
StatsAPI.params(x::ClusterPermutationTest) = x.cpc.stats
initial_fits(x::ClusterPermutationTest) = x.cpc.m
initial_fit(x::ClusterPermutationTest, i::Integer) = x.cpc.m[i]

use_threads_default(::ClusterPermutationTest) = true

function cluster_statistics(x::ClusterPermutationTest)
	return cluster_statistics(x.cpc.mass_fnc, x.cpc.stats, x.cpc.cc)
end;

function pvalues(x::ClusterPermutationTest;
	inhibit_warning::Bool = false)::Vector{Float64}
	# Monte Carlo permutation p value

	n = npermutations(x)
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
	spl = fits(x)
	clstats = cluster_statistics(x)
	for (i, stats) in enumerate(clstats)
		n_l = sum(spl[:, i] .> abs(stats))
		append!(p, 2 * (n_l / n))
	end
	return p
end;

function cluster_table(x::ClusterPermutationTest)
	cl_ranges = cluster_ranges(x)
	stats = params(x)
	pvals = pvalues(x; inhibit_warning = true)
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
		min = [minimum(stats[c]) for c in cl_ranges],
		max = [maximum(stats[c]) for c in cl_ranges],
		p = p,
		sign = sign)
end

function _test_info(x::ClusterPermutationTest)
	return "$(typeof(x)) ($(x.cpc.mass_fnc))"
end

function StatsAPI.summary(x::ClusterPermutationTest)
	println("$(_test_info(x))")
	println("  data: $(nepochs(x)) x $(epoch_length(x))")
	pt = pretty_table(String, cluster_table(x);
		show_subheader = false,
		formatters = (ft_printf("%0.2f", [5, 6]),
			ft_printf("%0.3f", [7])),
		vlines = :none,
		tf = tf_unicode_rounded)
	print(pt)
	return println("n permutations: $(npermutations(x))")
end;

function Base.show(io::IO, mime::MIME"text/plain", x::ClusterPermutationTest)
	clr = cluster_ranges(x)
	cc = cluster_criterium(x)
	println(io, "$(_test_info(x))")
	println(io, "  data: $(nepochs(x)) x $(epoch_length(x))")
	if cc isa ClusterDefinition
		println(io,
			"  $(length(clr)) cluster (ranges=$(cc.ranges)):")
	else
		println(io,
			"  $(length(clr)) cluster (threshold=$(cc.threshold), min_size=$(cc.min_size)):")
	end
	return println(io, "  $(npermutations(x)) permutations")
end;
