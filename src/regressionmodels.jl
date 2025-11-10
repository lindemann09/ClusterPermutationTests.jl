abstract type CPRegressionModel <: ClusterPermutationTest end

struct CPLinearModel <: CPRegressionModel
	cpc::CPCollection{StatsModels.TableRegressionModel}
	dat::CPData

	f::FormulaTerm
	effect::String # name of effect to be tested
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
end;

n_threads_default(::CPLinearModel) = Threads.nthreads()

####
#### Test info
####
function test_info(x::CPRegressionModel)
	return "$(typeof(x)) (effect=$(x.effect), $(x.cpc.mass_fnc))\n  $(x.f))"
end

####
#### Fit Regressions
####
function StatsAPI.fit(T::Type{<:CPRegressionModel},
	f::FormulaTerm, dat::CPData, cluster_criterium::TClusterCritODef; kwargs...)

	# use first rhs variable from formula as effect
	if f.rhs isa Term
		effect = f.rhs
	else
		effect = first(f.rhs)
	end
	fit(T, f, effect.sym, dat, cluster_criterium; kwargs...)
end

function StatsAPI.fit(::Type{<:CPLinearModel},
	f::FormulaTerm,
	effect::SymbolOString,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}())

	effect = String(effect)
	iv = first(split(effect, ": "))
	tbl = _prepare_design_table(f, dat.design, dv_dtype = eltype(dat.epochs))
	cpc = CPCollection{StatsModels.TableRegressionModel}(iv, mass_fnc, cluster_criterium)
	rtn = CPLinearModel(cpc,
		CPData(dat.epochs, tbl; unit_obs = unit_observation(dat.design)),
		f, effect, contrasts)
	fit_initial_time_series!(rtn)
	_check_effects(rtn.cpc.m[1], effect)
	return rtn
end


####
#### Sample statistics
####
function time_series_stats(cpt::CPRegressionModel, effect::SymbolOString)
	length(cpt.cpc.m) > 0 || return Float64[]
	# extract test statistics from initial fit
	i = _get_coefficient_row(cpt.cpc.m[1], String(effect))
	return time_series_stats(cpt::CPRegressionModel, i)
end

time_series_stats(cpt::CPRegressionModel, parameter_id::Integer) =
	get_parameter(cpt.cpc.ts, parameter_id)

####
#### Parameter estimates
####
"""estimates for a specific section in the time series (cluster) for a given permutation
"""
@inline function parameter_estimates(cpt::CPLinearModel,
	cl_ranges::Vector{TClusterRange};
	store_fits::Bool = false)::TParameterVector

	design = columntable(dat.design)
	param = TParameterVector()
	i = nothing # index for coefficient of iv
	dv_data = getproperty(design, cpt.f.lhs.sym)
	for (t, dv) in enumerate(eachcol(dat.epochs))
		dv_data[:] = dv
		md = fit(LinearModel, cpt.f, design; contrasts = cpt.contrasts) ## fit model!
		z = coef(md) ./ stderror(md) # parameter: z or t-value of effect
		if store_fits
			push!(cpt.cpc.m, md)
			push!(cpt.cpc.ts, TPStats(t, z))
		else
			#calc z for effect
			push!(param, TPStats(t, z))
		end
	end
	return param
end


##
## permutation stats
##
permutation_stats(cpt::CPRegressionModel, effect::Union{Integer, Symbol, String}) =
	_permutation_stats(cpt.cpc, cluster_ranges(cpt), effect)

##
## Cluster Functions
##
cluster_ranges(cpt::CPTTest, effect::Union{Integer, Symbol, String}) =
		_cluster_ranges(time_series_stats(cpt, effect), cpt.cpc.cc)

function cluster_mass(cpt::CPTTest, effect::Union{Integer, Symbol, String})
	ts = time_series_stats(cpt, effect)
	cl_ranges = _cluster_ranges(ts, cpt.cpc.cc)
	return _cluster_mass(cpt.cpc.mass_fnc, ts, cl_ranges)
end

function cluster_pvalues(cpt::CPTTest, effect::Union{Integer, Symbol, String};
		inhibit_warning::Bool = false)

	if effect isa Integer
		effect_id = effect
	else
		effect_id = _get_coefficient_row(cpt.cpc.m[1], String(effect))
	end
	return _cluster_pvalues(permutation_stats(cpt, effect_id),
				cluster_mass(cpt, effect_id), inhibit_warning)
end

cluster_table(cpt::CPTTest, effect::Union{Integer, Symbol, String}) =
	_cluster_table(time_series_stats(cpt), cluster_ranges(cpt, effect), cluster_pvalues(cpt, effect))


###
### Utilities for regression design tables
###

"""select columns from formula and add empty column for dependent variable"""
function _prepare_design_table(f::FormulaTerm, design::StudyDesign;
	dv_dtype::Type = Float64)::Table
	# get all predictors
	pred = Symbol[]
	_add_all_vars!(pred, f.rhs)

	for v in pred
		has_variable(design, v) || throw(ArgumentError("Variable '$(v)' not found in design table!"))
	end

	perm_design = select_col(columntable(design), pred)  # select required variables
	dv = Vector{dv_dtype}(undef, length(design))
	dv_name = Symbol(f.lhs.sym)
	return Table(perm_design, (; dv_name => dv)) # add dependent variable column
end

function _check_effects(model::RegressionModel, effect::String)
	n = coefnames(model)
	try
		i = _get_coefficient_row(model, effect)
		if i == 0
			@warn("Effect '$effect' not found in model. Effects: '$(n)'.")
		end
	catch e
		@warn("Effect '$effect' unclear. Effects: '$(n)'.") # found multiple effects
	end
end

function _get_coefficient_row(md::StatisticalModel, effect::String)
	# get row id of a specific effect in the coefficient vector/coeftable
	# for binary  categorial variables is the variable names is sufficient and level (..: level) is not required

	# find exact match
	n = coefnames(md)
	i = findfirst(x -> x == effect, n)
	if isnothing(i)
		# find binary categorical variable match
		idx = findall(x -> startswith(x, "$effect: "), n)
		if length(idx) == 1
			i = idx[1] # only two levels
		elseif length(idx) > 1
			throw(ArgumentError("Effect '$effect' unclear. Available effects: '$(n)'."))
		else
			return 0
		end
	end
	i > 0 || throw(ArgumentError("Can not find effect '$effect'."))
	return i
end

function _add_all_vars!(vec::Vector{Symbol}, x::Tuple)
	for t in x
		_add_all_vars!(vec, t)
	end
	return vec
end

_add_all_vars!(vec::Vector{Symbol}, ::ConstantTerm) = vec # do nothing
_add_all_vars!(vec::Vector{Symbol}, t::Term) = t.sym in vec ? vec : push!(vec, t.sym)
_add_all_vars!(vec::Vector{Symbol}, x::InteractionTerm) = _add_all_vars!(vec, x.terms)
_add_all_vars!(vec::Vector{Symbol}, x::FunctionTerm) = _add_all_vars!(vec, Tuple(x.args))
