###
### AbstractClusterPermutationTest
###

const no_effect_error = ArgumentError("The model has multiple coefficients. Please specify an effect.")

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


"""
    time_series_fits(x::ClusterPermutationTest)

Return the vector of fitted models from the initial (un-permuted) fit, one per time point.
"""
time_series_fits(x::ClusterPermutationTest) = x.cpc.Md

"""
    npermutations(x::ClusterPermutationTest)

Return the number of permutations accumulated so far (via `resample!`). Returns 0 if
`resample!` has not yet been called.
"""
function npermutations(x::ClusterPermutationTest)
	return length(x.cpc.X)
end

"""
	ncoefs(x::ClusterPermutationTest)

Return the number of coefficients in the model, i.e., the number of effects
"""
ncoefs(x::ClusterPermutationTest) = size(x.cpc.coefs, 2)

"""
    time_series_stats(x::ClusterPermutationTest, [effect])

Return the time series of test statistics for the specified `effect` from the initial fit.

`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name
(see `coefnames`).
"""
function time_series_stats(x::ClusterPermutationTest)
	ncoefs(x) == 1 ? time_series_stats(x, 1) : throw(no_effect_error)
end
time_series_stats(x::ClusterPermutationTest, effect::Union{Integer, Symbol, String}) = view(x.cpc.coefs, :, _effect_id(x, effect))

##
## Cluster Functions
##


"""
    cluster(cpt::ClusterPermutationTest, effect)

Return the detected or defined cluster ranges for the specified `effect`.

`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
function cluster(x::ClusterPermutationTest)
	ncoefs(x) == 1 ? cluster(x, 1) : throw(no_effect_error)
end
function cluster(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String})
	ts = view(cpt.cpc.coefs, :, _effect_id(cpt, effect)) # time series stats for this effect
	return _cluster_ranges(ts, cluster_type(cpt.config))
end

"""
    cluster_mass_stats(cpt::ClusterPermutationTest, effect)

Return the cluster mass statistics for each detected or defined cluster of the specified `effect`.

The mass statistic is computed by applying `mass_fnc` (default: `sum`) to the
time-series statistics within each cluster range.
`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
function cluster_mass_stats(x::ClusterPermutationTest)
	ncoefs(x) == 1 ? cluster_mass_stats(x, 1) : throw(no_effect_error)
end
function cluster_mass_stats(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String})
	i = _effect_id(cpt, effect)
	ts = time_series_stats(cpt, i)
	cl_ranges = cluster(cpt, i)
	return _cluster_mass_stats(cpt.config.mass_fnc, ts, cl_ranges)
end

"""
    cluster_pvalues(cpt::ClusterPermutationTest, effect; inhibit_warning=false, clusterwise=false)

Return the Monte Carlo permutation p-values for each detected cluster of the specified `effect`.

Requires at least 1000 permutations; warns when fewer than 5000 are available.
Call `resample!` to accumulate permutations.
`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
function cluster_pvalues(x::ClusterPermutationTest; kwargs...)
	ncoefs(x) == 1 ? cluster_pvalues(x, 1; kwargs...) : throw(no_effect_error)
end
function cluster_pvalues(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String};
	inhibit_warning::Bool = false,
	clusterwise::Bool = false)
	i = _effect_id(cpt, effect)
	return _cluster_pvalues(cluster_nhd(cpt, i; clusterwise), cluster_mass_stats(cpt, i), inhibit_warning)
end

"""
    cluster_table(cpt::ClusterPermutationTest; add_effect_names=false, one_tail=false, clusterwise=false)
    cluster_table(cpt::ClusterPermutationTest, effect; inhibit_warning=false, add_effect_names=false,
				one_tail=false, clusterwise=false)

Return a table summarising detected clusters with their range, size, mass statistic, and p-value.

When called without an `effect`, results for all effects are combined into a single table.
`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
function cluster_table(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String};
	inhibit_warning::Bool = false,
	add_effect_names::Bool = false,
	one_tail::Bool = false,
	clusterwise::Bool = false)::CoefTable

	i = _effect_id(cpt, effect)
	coef_name = coefnames(cpt)[i]
	ts = time_series_stats(cpt, i)
	cl_ranges = cluster(cpt, i)
	cl_mass_stats = _cluster_mass_stats(cpt.config.mass_fnc, ts, cl_ranges)
	p_vals = _cluster_pvalues(cluster_nhd(cpt, i; clusterwise), cl_mass_stats, inhibit_warning; one_tail)
	return _cluster_table(i, coef_name, cl_ranges, cl_mass_stats, p_vals; add_effect_names)
end

"""Cluster table for all effects"""
function cluster_table(cpt::ClusterPermutationTest; kwargs...)::CoefTable

	rtn = cluster_table(cpt, 1; kwargs...)
	# add all other effects
	for eid in 2:ncoefs(cpt)
		tmp = cluster_table(cpt, eid; one_tail, add_effect_names, inhibit_warning = true)
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
    cluster_nhd(cpt::ClusterPermutationTest, [effect]; clusterwise::Bool=false)

Return the bootstrapped null-hypothesis distribution (NHD; i.e., permutation distribution) of
cluster mass statistics for the specified `effect`. The function returns the NHD based on the
permutation distribution defined in `cpt.config.permutation_distr` (e.g. maxmass, cluster-mass).

Each column is the NHD of for cluster statistic of the detected cluster in this effect. If `maxmass`
permutation distribution was used, columns are identical, since all clusters are tested against
the same NHD.

If `clusterwise` if true, the cluster-wise NHDs are return, even if `massmass` was used while
model fitting.

The returned matrix has shape `(n_permutations × n_clusters)`. Returns an empty matrix if
`resample!` has not been called yet.
`effect` can be an integer index, a `Symbol`, or a `String` matching a coefficient name.
"""
function cluster_nhd(x::ClusterPermutationTest; kwargs...)
	ncoefs(x) == 1 ? cluster_nhd(x, 1; kwargs...) : throw(no_effect_error)
end
function cluster_nhd(cpt::ClusterPermutationTest,
	effect::Union{Integer, Symbol, String};
	clusterwise::Bool=false)::TParameterMatrix # (permutation X cluster)

	if length(cpt.cpc.X) == 0
		return zeros(Float64, 0, 0)
	else
		e_id = _effect_id(cpt, effect)
		n_cluster = length(cluster(cpt, e_id))
		if is_maxmass(cpt.config) && !clusterwise
			mm = [x.max_mass[e_id] for x in cpt.cpc.X]
			return mm * ones(Float64, 1, n_cluster) # repeat the same column for all clusters
		else
			return reduce(hcat, [x.cluster_mass[e_id] for x in cpt.cpc.X])'
		end
	end
end

function Base.summary(x::ClusterPermutationTest)
	println(_info(x))
	ivs = join(string.(x.cpc.shuffle_ivs), ", ")
	println("  shuffled variables: $(ivs)")
	println("  cluster $(_cluster_info_str(x.config.cc))")
	println("  cluster stats: $(x.config.mass_fnc), statistic: $(x.config.permutation_distr)")
	display(cluster_table(x))
	return println("  n permutations: $(npermutations(x))")
end;

function Base.show(io::IO, mime::MIME"text/plain", x::ClusterPermutationTest)
	println(io, _info(x))
	return println(io, "  $(npermutations(x)) permutations")
end;

function _info(x::ClusterPermutationTest)::String
	rtn = "$(test_info(x))\n"
	rtn *= "  data: $(nepochs(x)) x $(epoch_length(x))\n"
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