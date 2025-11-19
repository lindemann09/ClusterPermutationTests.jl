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

	# shuffle_ivs: no covariates and no random effects (only categorical predictors)
	pred = StatsModels.termvars(f.rhs) # get all predictor variables
	i = findall(is_randomeffectsterm.(f.rhs))
	random_effects = isempty(i) ? Symbol[] : [x.args[2].sym for x in f.rhs[i]]
	shuffle_ivs = filter(x -> !is_covariate(dat.design, x) && x ∉ random_effects, pred)

	fit(T, f, shuffle_ivs, dat, cluster_criterium; kwargs...)
end

function StatsAPI.fit(::Type{<:CPLinearModel},
	f::FormulaTerm,
	shuffle_ivs::Union{Vector{Symbol}, Symbol, Vector{String}, String},
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}())

	data, shuffle_ivs = _prepare_regression_data(f, dat, shuffle_ivs)
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
	design::AbstractStudyDesign,
	time_points::Vector{Int32};
	store_fits::Bool = false)::T2DParamVector # time x effect

	design = columntable(design)
	param = T2DParamVector()
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

function _prepare_regression_data(f::FormulaTerm, dat::CPData,
	shuffle_ivs::Union{Vector{Symbol}, Symbol, Vector{String}, String})::Tuple{CPData, Vector{Symbol}}

	tbl = __prepare_design_table(f, dat.design; dv_dtype = eltype(dat.epochs))
	data = CPData(dat.epochs, tbl; unit_obs = unit_observation(dat.design))

	shuffle_ivs = StudyDesigns.to_symbol_vector(shuffle_ivs)
	for v in shuffle_ivs
		has_variable(data.design, v) || throw(ArgumentError("Variable '$(v)' not found in formula or design table!"))
	end
	return data, shuffle_ivs
end

"""select columns from formula and add (optionally) empty column for dependent variable"""
function __prepare_design_table(f::FormulaTerm, design::AbstractStudyDesign;
	dv_dtype::Type = Float64)::Table

	pred = StatsModels.termvars(f.rhs) # get all predictor variables
	for v in pred
		has_variable(design, v) || throw(ArgumentError("Variable '$(v)' not found in design table!"))
	end

	perm_design = select_cols(columntable(design), pred)  # select required variables

	# add_dv_column
	dv_name = f.lhs.sym
	dv = zeros(dv_dtype, length(design))
	return Table(perm_design, (; dv_name => dv)) # add dependent variable column
end