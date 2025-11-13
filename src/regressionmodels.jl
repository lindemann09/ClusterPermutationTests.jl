abstract type CPRegressionModel <: ClusterPermutationTest end

struct CPLinearModel <: CPRegressionModel
	cpc::CPCollection{StatsModels.TableRegressionModel}
	dat::CPData

	f::FormulaTerm
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
end;

n_threads_default(::CPLinearModel) = Threads.nthreads()

####
#### Test info
####
function test_info(x::CPRegressionModel)
	return "$(typeof(x)) ($(x.cpc.mass_fnc))\n  $(x.f)"
end

####
#### Fit Regressions
####
function StatsAPI.fit(T::Type{<:CPRegressionModel},
	f::FormulaTerm, dat::CPData, cluster_criterium::TClusterCritODef; kwargs...)
	# default shuffle variables: all categorical predictors except covariates and random effects

	shuffle_ivs = _all_predictors(f; skip_randeff = true)
	shuffle_ivs = filter(x->!is_covariate(dat.design, x), shuffle_ivs) # no covariates (only categorical)
	fit(T, f, shuffle_ivs, dat, cluster_criterium; kwargs...)
end

function StatsAPI.fit(::Type{<:CPLinearModel},
	f::FormulaTerm,
	shuffle_ivs::Union{Vector{Symbol}, Symbol, Vector{String}, String},
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}())

	tbl = _prepare_design_table(f, dat.design, dv_dtype = eltype(dat.epochs))
	data = CPData(dat.epochs, tbl; unit_obs = unit_observation(dat.design))

	shuffle_ivs = StudyDesigns.to_symbol_vector(shuffle_ivs)
	for v in shuffle_ivs
		has_variable(data.design, v) || throw(ArgumentError("Variable '$(v)' not found in formula or design table!"))
	end

	cpc = CPCollection{StatsModels.TableRegressionModel}(shuffle_ivs, mass_fnc, cluster_criterium)
	rtn = CPLinearModel(cpc, data, f, contrasts)
	fit_initial_time_series!(rtn)
	return rtn
end

####
#### coefnames
####
function StatsAPI.coefnames(cpt::CPRegressionModel)
	rtn = coefnames(first(cpt.cpc.m))
	return rtn[2:end] # remove Intercept
end

####
#### Parameter estimates
####
"""estimates for a specific section in the time series (cluster) for a given permutation
"""
@inline function parameter_estimates(cpt::CPLinearModel,
	design::AbstractStudyDesign;
	fit_cluster_only::Bool = true,
	store_fits::Bool = false)::TVecTimeXParameter

	design = columntable(design)
	param = TVecTimeXParameter()
	if fit_cluster_only
		time_points = cpt.cpc.tp
	else
		time_points = Int32(1):Int32(epoch_length(cpt.dat))
	end

	dv_data = getproperty(design, cpt.f.lhs.sym)
	for t in time_points
		dv_data[:] = cpt.dat.epochs[:, t] # update dependent variable FIXME view?
		md = fit(LinearModel, cpt.f, design; contrasts = cpt.contrasts) ## fit model!
		z = coef(md) ./ stderror(md) # parameter: z or t-value of effect
		push!(param, z[2:end])
		if store_fits
			push!(cpt.cpc.m, md)
		end
	end
	return param
end


###
### Utilities for regression design tables
###

"""select columns from formula and add empty column for dependent variable"""
function _prepare_design_table(f::FormulaTerm, design::AbstractStudyDesign;
	dv_dtype::Type = Float64)::Table
	# get all predictors
	pred = _all_predictors(f)
	for v in pred
		has_variable(design, v) || throw(ArgumentError("Variable '$(v)' not found in design table!"))
	end

	perm_design = select_col(columntable(design), pred)  # select required variables
	dv = Vector{dv_dtype}(undef, length(design))
	dv_name = Symbol(f.lhs.sym)
	return Table(perm_design, (; dv_name => dv)) # add dependent variable column
end

function _all_predictors(f::FormulaTerm; skip_randeff::Bool = false)::Vector{Symbol}
	pred = Symbol[]
	_add_all_vars!(pred, f.rhs, skip_randeff)
end

function _add_all_vars!(vec::Vector{Symbol}, x::Tuple, skip_randeff)
	for t in x
		_add_all_vars!(vec, t, skip_randeff)
	end
	return vec
end

_add_all_vars!(vec::Vector{Symbol}, ::ConstantTerm, ::Bool) = vec # do nothing
_add_all_vars!(vec::Vector{Symbol}, t::Term, ::Bool) = t.sym in vec ? vec : push!(vec, t.sym)
_add_all_vars!(vec::Vector{Symbol}, x::InteractionTerm, sr::Bool) = _add_all_vars!(vec, x.terms, sr)
function _add_all_vars!(vec::Vector{Symbol}, x::FunctionTerm, skip_randeff::Bool)
	if !skip_randeff || !is_randomeffectsterm(x)
		_add_all_vars!(vec, Tuple(x.args), skip_randeff)
	end
end