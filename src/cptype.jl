###
### CPCollection
###
const TParameterVector = Vector{Float64}
const TParameterMatrix = Matrix{Float64}
const T2DParamVector = Vector{TParameterVector}
const no_effect_error = ArgumentError("Please specify an effect.")

mutable struct CPCollection{M}
	shuffle_ivs::Vector{Symbol} # name of the to be shuffled independent variable
	mass_fnc::Function # cluster mass function
	cc::TClusterCritODef # cluster definition

	M::Vector{M} # fitted models of initial fit
	coefs::TParameterMatrix # (time X effect) time series statistics of the initial fit

	cl::Vector{Vector{TClusterRange}} # cluster ranges for all effects (effectXcluster)

	S::Vector{TParameterMatrix} # one matrix per effect, (permutation X cluster)
end;

function CPCollection{M}(shuffle_ivs::Vector{Symbol}, mass_fnc::Function,
	cluster_criterium::TClusterCritODef) where {M}
	return CPCollection{M}(shuffle_ivs, mass_fnc, cluster_criterium,
			M[], zeros(Float64, 0, 0), TClusterRange[], TParameterMatrix[])
end

###
### ClusterPermutationTest
###

"""
Abstract base type for all cluster permutation tests.

Concrete subtypes (e.g. `CPPairedSampleTTest`, `CPLinearModel`) are constructed via
`fit(T, ...)`, which runs the initial time-series fit and detects clusters. Call
`resample!` afterwards to build the null-hypothesis distribution.
"""
abstract type ClusterPermutationTest end
#requires
#	cpc::CPCollection
#	dat::CPData

nepochs(x::ClusterPermutationTest) = nepochs(x.dat)
epoch_length(x::ClusterPermutationTest) = epoch_length(x.dat)
design_table(x::ClusterPermutationTest) = design_table(x.dat)
StudyDesigns.unit_observation(x::ClusterPermutationTest) = unit_observation(x.dat.design.uo)

cluster_criterium(x::ClusterPermutationTest) = x.cpc.cc

"""
    time_series_fits(x::ClusterPermutationTest)

Return the vector of fitted models from the initial (un-permuted) fit, one per time point.
"""
time_series_fits(x::ClusterPermutationTest) = x.cpc.M

"""
    npermutations(x::ClusterPermutationTest)

Return the number of permutations accumulated so far (via `resample!`). Returns 0 if
`resample!` has not yet been called.
"""
function npermutations(x::ClusterPermutationTest)
	if length(x.cpc.S) == 0
		return 0
	else
		return size(x.cpc.S[1], 1)
	end
end
ncoefs(x::ClusterPermutationTest) = size(x.cpc.coefs, 2)

"""
    time_series_stats(x::ClusterPermutationTest, effect)

Return the time series of test statistics for the specified `effect` from the initial fit.

`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name
(see `coefnames`).
"""
time_series_stats(::ClusterPermutationTest) = throw(no_effect_error)
time_series_stats(x::ClusterPermutationTest, effect::Union{Integer, Symbol, String}) = view(x.cpc.coefs, :, _effect_id(x, effect))

##
## Cluster Functions
##

"""
    cluster(cpt::ClusterPermutationTest, effect)

Return the detected cluster ranges for the specified `effect`.

`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
cluster(::ClusterPermutationTest) = throw(no_effect_error)
cluster(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String}) = cpt.cpc.cl[_effect_id(cpt, effect)]

"""
    cluster_mass_stats(cpt::ClusterPermutationTest, effect)

Return the cluster mass statistics for each detected cluster of the specified `effect`.

The mass statistic is computed by applying `mass_fnc` (default: `sum`) to the
time-series statistics within each cluster range.
`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
cluster_mass_stats(::ClusterPermutationTest) = throw(no_effect_error)
function cluster_mass_stats(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String})
	i = _effect_id(cpt, effect)
	ts = time_series_stats(cpt, i)
	cl_ranges = cluster(cpt, i)
	return _cluster_mass_stats(cpt.cpc.mass_fnc, ts, cl_ranges)
end

"""
    cluster_pvalues(cpt::ClusterPermutationTest, effect; inhibit_warning=false)

Return the Monte Carlo permutation p-values for each detected cluster of the specified `effect`.

Requires at least 1000 permutations; warns when fewer than 5000 are available.
Call `resample!` to accumulate permutations.
`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
cluster_pvalues(::ClusterPermutationTest; kwargs...) = throw(no_effect_error)
function cluster_pvalues(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String};
	inhibit_warning::Bool = false)
	i = _effect_id(cpt, effect)
	return _cluster_pvalues(cluster_nhd(cpt, i), cluster_mass_stats(cpt, i), inhibit_warning)
end

"""
    cluster_table(cpt::ClusterPermutationTest)
    cluster_table(cpt::ClusterPermutationTest, effect; inhibit_warning=false, add_effect_names=false)

Return a table summarising detected clusters with their range, size, mass statistic, and p-value.

When called without an `effect`, results for all effects are combined into a single table.
`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
function cluster_table(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String};
	inhibit_warning::Bool = false,
	add_effect_names::Bool = false)::CoefTable
	i = _effect_id(cpt, effect)
	coef_name = coefnames(cpt)[i]
	ts = time_series_stats(cpt, i)
	cl_ranges = cluster(cpt, i)
	cl_mass_stats = _cluster_mass_stats(cpt.cpc.mass_fnc, ts, cl_ranges)
	p_vals = _cluster_pvalues(cluster_nhd(cpt, i), cl_mass_stats, inhibit_warning)
	return _cluster_table(i, coef_name, cl_ranges, cl_mass_stats, p_vals; add_effect_names)
end

"""Cluster table for all effects"""
function cluster_table(cpt::ClusterPermutationTest)::CoefTable
	add_effect_names = ncoefs(cpt) > 1
	rtn = cluster_table(cpt, 1; add_effect_names)
	# add all other effects
	for eid in 2:ncoefs(cpt)
		tmp = cluster_table(cpt, eid; add_effect_names, inhibit_warning = true)
		for (x, y) in zip(rtn.cols, tmp.cols)
			append!(x, y)
		end
		append!(rtn.rownms, tmp.rownms)
	end
	return rtn
end

##
## Null-hypothesis distributions
##
"""
    cluster_nhd(cpt::ClusterPermutationTest, effect)

Return the null-hypothesis distribution of cluster mass statistics for the specified `effect`.

The returned matrix has shape `(n_permutations × n_clusters)`. Returns an empty matrix if
`resample!` has not been called yet.
`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
cluster_nhd(::ClusterPermutationTest) = throw(no_effect_error)
function cluster_nhd(cpt::ClusterPermutationTest,
	effect::Union{Integer, Symbol, String})::TParameterMatrix # (permutation X cluster)
	if length(cpt.cpc.S) == 0
		return zeros(Float64, 0, 0)
	else
		e_id = _effect_id(cpt, effect)
		return cpt.cpc.S[e_id]
	end
end

function Base.summary( x::ClusterPermutationTest)
	println(_info(x))
	display(cluster_table(x))
	return println("  n permutations: $(npermutations(x))")
end;

function Base.show(io::IO, mime::MIME"text/plain", x::ClusterPermutationTest)
	println(io, _info(x))
	cc = cluster_criterium(x)
	if cc isa ClusterDefinition
		println(io, "  cluster definition: ranges=$(cc.ranges)")
	else
		println(io,
			"  cluster definition: threshold=$(cc.threshold), min_size=$(cc.min_size)")
	end
	return println(io, "  $(npermutations(x)) permutations")
end;

function _info(x::ClusterPermutationTest)::String
	rtn = "$(test_info(x))\n"
	rtn *= "  data: $(nepochs(x)) x $(epoch_length(x))\n"
	ivs = join(string.(x.cpc.shuffle_ivs), ", ")
	rtn *= "  shuffled ivs: $(ivs)\n"
	n = join(coefnames(x), "\n           ")
	return rtn * "  effects: $n"
end

####
#### Helper functions
####
function _effect_id(cpt::ClusterPermutationTest, effect::Integer)
	(effect > ncoefs(cpt) || effect < 1) &&
		throw(ArgumentError("Effect index $(effect) out of bounds."))
	return effect
end
function _effect_id(cpt::ClusterPermutationTest, effect::Union{Symbol, String})
	names = coefnames(cpt)
	rtn = findfirst(isequal(string(effect)), names)
	rtn === nothing &&
		throw(ArgumentError("Effect '$(effect)' not found in model coefficients. Used names: '$(names)'."))
	return rtn
end