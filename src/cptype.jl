###
### CPCollection
###
const TParameterVector = Vector{Float64}
const TParameterMatrix = Matrix{Float64}
const TVecTimeXParameter = Vector{TParameterVector}
const no_effect_error = ArgumentError("Please specify an effect.")

mutable struct CPCollection{M}
	iv::Symbol # name of the to be shuffled independent variable
	mass_fnc::Function # cluster mass function
	cc::TClusterCritODef # cluster definition

	m::Vector{M} # fitted models of initial fit
	coefs::TParameterMatrix # (time X effect) time series statistics of the initial fit

	tp::Vector{Int32}# time points involved in all clusters
	S::Vector{TParameterMatrix} # one matrix per effect, (time X permutation)
end;

CPCollection{M}(iv::SymbolOString, mass_fnc::Function, cluster_criterium::TClusterCritODef) where {M} =
	CPCollection(Symbol(iv), mass_fnc, cluster_criterium,
		M[], zeros(Float64, 0, 0), Int32[], TParameterMatrix[])

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
time_series_stats(x::ClusterPermutationTest, effect::Union{Integer, Symbol, String}) =
	view(x.cpc.coefs, :, _effect_id(x, effect))
time_series_stats(::ClusterPermutationTest) = throw(no_effect_error)

ncoefs(x::ClusterPermutationTest) = size(x.cpc.coefs, 2)

##
## permutation stats
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
		for cl in cl_ranges
			rows = findfirst(isequal(cl.start), cpt.cpc.tp):findfirst(isequal(cl.stop), cpt.cpc.tp)
			mtx_sub = view(cpt.cpc.S[e_id], rows, :)
			cl_mass = [cpt.cpc.mass_fnc(col) for col in eachcol(mtx_sub)]
			push!(rtn, cl_mass)
		end
		return reduce(hcat, rtn)
	end
end

##
## Cluster Functions
##

cluster_ranges(::ClusterPermutationTest) = throw(no_effect_error)
cluster_ranges(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String}) =
	_cluster_ranges(time_series_stats(cpt, effect), cpt.cpc.cc) #FIXME effect strings

cluster_mass_stats(::ClusterPermutationTest) = throw(no_effect_error)
function cluster_mass_stats(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String})
	ts = time_series_stats(cpt, effect)
	cl_ranges = _cluster_ranges(ts, cpt.cpc.cc)
	return _cluster_mass_stats(cpt.cpc.mass_fnc, ts, cl_ranges)
end

cluster_pvalues(::ClusterPermutationTest; kwargs...) = throw(no_effect_error)
function cluster_pvalues(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String};
	inhibit_warning::Bool = false)
	i = _effect_id(cpt, effect)
	return _cluster_pvalues(cluster_nhd(cpt, i), cluster_mass_stats(cpt, i), inhibit_warning)
end

function cluster_table(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String};
	inhibit_warning::Bool = false)::CoefTable
	i = _effect_id(cpt, effect)
	coef_name = coefnames(cpt)[i]
	ts = time_series_stats(cpt, i)
	cl_ranges = _cluster_ranges(ts, cpt.cpc.cc)
	cl_mass_stats = _cluster_mass_stats(cpt.cpc.mass_fnc, ts, cl_ranges)
	p_vals = _cluster_pvalues(cluster_nhd(cpt, i), cl_mass_stats, inhibit_warning)
	_cluster_table(i, coef_name, cl_ranges, cl_mass_stats, p_vals)
end

function cluster_table(cpt::ClusterPermutationTest)::CoefTable
	rtn = cluster_table(cpt, 1)
	# add all other effects
	for eid in 2:ncoefs(cpt)
		tmp = cluster_table(cpt, eid; inhibit_warning = true)
		for (x, y) in zip(rtn.cols, tmp.cols)
       		append!(x, y)
       	end
		append!(rtn.rownms, tmp.rownms)
	end
	return rtn
end


function _effect_id(cpt::ClusterPermutationTest, effect::Union{Integer, Symbol, String})
	if effect isa Integer
		(effect > ncoefs(cpt) || effect < 1) &&
			throw(ArgumentError("Effect index $(effect) out of bounds."))
		return effect
	else
		names = coefnames(cpt)
		rtn = findfirst(isequal(String(effect)), names)
		rtn === nothing &&
			throw(ArgumentError("Effect '$(effect)' not found in model coefficients. Used names: '$(names)'."))
		return rtn
	end
end

StatsAPI.summary(x::ClusterPermutationTest) = summary(x, 1)
function StatsAPI.summary(x::ClusterPermutationTest, effect::Union{Integer, Symbol, String})
	println("$(test_info(x))")
	println("  data: $(nepochs(x)) x $(epoch_length(x))")
	eid = _effect_id(x, effect)
	println("  Effect '$(coefnames(x)[eid])'")
	pt = pretty_table(String, cluster_table(x, eid);
		show_subheader = false,
		formatters = (ft_printf("%0.2f", [5, 6]),
			ft_printf("%0.3f", [7])),
		vlines = :none,
		tf = tf_unicode_rounded)
	print(pt)
	return println("n permutations: $(npermutations(x))")
end;

function Base.show(io::IO, mime::MIME"text/plain", x::ClusterPermutationTest)
	println(io, "$(test_info(x))")
	println(io, "  data: $(nepochs(x)) x $(epoch_length(x))")
	clr = TClusterRange[]
	cc = cluster_criterium(x)
	if cc isa ClusterDefinition
		println(io, "  cluster definition: ranges=$(cc.ranges)")
	else
		println(io,
			"  cluster definition: threshold=$(cc.threshold), min_size=$(cc.min_size)")
	end
	return println(io, "  $(npermutations(x)) permutations")
end;
