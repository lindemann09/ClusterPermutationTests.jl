###
### CPCollection
###
const TParameterVector = Vector{Float64}
const TParameterMatrix = Matrix{Float64}
const TVecTimeXParameter = Vector{TParameterVector}
const no_effect_error = ArgumentError("Please specify an effect.")

mutable struct CPCollection{M}
	shuffle_ivs::Vector{Symbol} # name of the to be shuffled independent variable
	mass_fnc::Function # cluster mass function
	cc::TClusterCritODef # cluster definition

	m::Vector{M} # fitted models of initial fit
	coefs::TParameterMatrix # (time X effect) time series statistics of the initial fit

	cl::Vector{Vector{TClusterRange}} # cluster ranges for all effects

	S::Vector{TParameterMatrix} # one matrix per effect, (time X permutation)
end;

CPCollection{M}(shuffle_ivs::Vector{Symbol}, mass_fnc::Function, cluster_criterium::TClusterCritODef) where {M} =
	CPCollection{M}(shuffle_ivs, mass_fnc, cluster_criterium,
		M[], zeros(Float64, 0, 0), TClusterRange[], TParameterMatrix[])

###
### ClusterPermutationTest
###

abstract type ClusterPermutationTest end
#requires
#	cpc::CPCollection
#	dat::CPData

nepochs(x::ClusterPermutationTest) = nepochs(x.dat)
epoch_length(x::ClusterPermutationTest) = epoch_length(x.dat)
design_table(x::ClusterPermutationTest) = design_table(x.dat)
StudyDesigns.unit_observation(x::ClusterPermutationTest) = unit_observation(x.dat.design.uo)

cluster_criterium(x::ClusterPermutationTest) = x.cpc.cc
time_series_fits(x::ClusterPermutationTest) = x.cpc.m

function npermutations(x::ClusterPermutationTest)
	if length(x.cpc.S) == 0
		return 0
	else
		return size(x.cpc.S[1], 2)
	end
end
ncoefs(x::ClusterPermutationTest) = size(x.cpc.coefs, 2)

time_series_stats(::ClusterPermutationTest) = throw(no_effect_error)
time_series_stats(x::ClusterPermutationTest, effect::Union{Integer, Symbol, String}) =
	view(x.cpc.coefs, :, _effect_id(x, effect))

##
## Cluster Functions
##

cluster_ranges(::ClusterPermutationTest) = throw(no_effect_error)
cluster_ranges(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String}) =
	cpt.cpc.cl[_effect_id(cpt, effect)]

cluster_mass_stats(::ClusterPermutationTest) = throw(no_effect_error)
function cluster_mass_stats(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String})
	i = _effect_id(cpt, effect)
	ts = time_series_stats(cpt, i)
	cl_ranges = cluster_ranges(cpt, i)
	return _cluster_mass_stats(cpt.cpc.mass_fnc, ts, cl_ranges)
end

cluster_pvalues(::ClusterPermutationTest; kwargs...) = throw(no_effect_error)
function cluster_pvalues(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String};
	inhibit_warning::Bool = false)
	i = _effect_id(cpt, effect)
	return _cluster_pvalues(cluster_nhd(cpt, i), cluster_mass_stats(cpt, i), inhibit_warning)
end

function cluster_table(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String};
	inhibit_warning::Bool = false,
	add_effect_names::Bool = false)::CoefTable

	i = _effect_id(cpt, effect)
	coef_name = coefnames(cpt)[i]
	ts = time_series_stats(cpt, i)
	cl_ranges = cluster_ranges(cpt, i)
	cl_mass_stats = _cluster_mass_stats(cpt.cpc.mass_fnc, ts, cl_ranges)
	p_vals = _cluster_pvalues(cluster_nhd(cpt, i), cl_mass_stats, inhibit_warning)
	_cluster_table(i, coef_name, cl_ranges, cl_mass_stats, p_vals; add_effect_names)
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
"""Null-hypothesis distributions of the cluster mass statistics"""
cluster_nhd(::ClusterPermutationTest) = throw(no_effect_error)
function cluster_nhd(cpt::ClusterPermutationTest,
	effect::Union{Integer, Symbol, String})::TParameterMatrix # (permutation X cluster)

	if length(cpt.cpc.S) == 0
		return zeros(Float64, 0, 0)
	else
		rtn = Vector{Float64}[]
		e_id = _effect_id(cpt, effect)
		cl_ranges = cluster_ranges(cpt, e_id)
		time_points = _joined_ranges(cpt.cpc.cl)
		for cl in cl_ranges
			rows = findfirst(isequal(cl.start), time_points):findfirst(isequal(cl.stop), time_points)
			mtx_sub = view(cpt.cpc.S[e_id], rows, :)
			cl_mass = [cpt.cpc.mass_fnc(col) for col in eachcol(mtx_sub)]
			push!(rtn, cl_mass)
		end
		return reduce(hcat, rtn)
	end
end


function StatsAPI.summary(x::ClusterPermutationTest)
	println(_info(x))
	println(cluster_table(x))
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