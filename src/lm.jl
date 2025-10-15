abstract type CPRegressionModel <: ClusterPermutationTest end

struct CPLinearModel <: CPRegressionModel
	cpc::CPCollection{StatsModels.TableRegressionModel}
	dat::CPData

	f::FormulaTerm
	effect::String # name of effect to be tested
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
end;

function StatsAPI.fit(T::Type{<:CPRegressionModel},
	f::FormulaTerm, dat::CPData, cluster_criterium::TClusterCritODef; kwargs...)

	# use first rhs variable from formula as effect
	if f.rhs isa Term
    	effect = f.rhs
	else
		effect = first(f.rhs)
		@warn("Effect is not defined. Using '$effect'.")
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


	if is_mixedmodel(f)
		throw(ArgumentError("Use CPMixedModels for mixed models."))
	end

	effect = String(effect)
	iv = first(split(effect, ": "))
	tbl = _prepare_design_table(f, dat.design, dv_dtype = eltype(dat.mtx))

	unit_obs = unit_observation(dat.design)
	if unit_obs ∉ columnnames(tbl)
		unit_obs = nothing
	end

	cpc = CPCollection{StatsModels.TableRegressionModel}(iv, mass_fnc, cluster_criterium)
	rtn = CPLinearModel(cpc, CPData(dat.mtx, tbl; unit_obs), f, effect, contrasts)
	initial_fit!(rtn)
	return rtn
end

####
#### definition of parameter_estimates
####

"""estimates for a specific section in the time series (cluster) for a given permutation
"""
@inline function parameter_estimates(cpt::CPLinearModel, dat::CPData;
		initial_fit::Bool = false)::TParameterVector
	dv_name = Symbol(cpt.f.lhs.sym)
	design = columntable(dat.design)
	param = TParameterVector() # TODO would be vector preallocation faster?
	i = nothing # index for coefficient of iv
	dv_data = getproperty(design, dv_name)
	for dv in eachcol(dat.mtx)
		dv_data[:] = dv
		md = fit(LinearModel, cpt.f, design; contrasts = cpt.contrasts) ## fit model!
		if isnothing(i)
			i = get_coefficient_row(md, cpt.effect)
		end
		z = coef(md)[i] / stderror(md)[i] # parameter: t-value of effect
		if initial_fit
			push!(cpt.cpc.m, md)
			push!(cpt.cpc.stats, z)
		else
			push!(param, z)
		end
	end
	return param
end

###
### Utilities for regression design tables
###

"""select columns from formula and add empty column for dependent variable"""
function _prepare_design_table(f::FormulaTerm, design::StudyDesign;
	 dv_dtype::Type=Float64)::Table

	dv = Vector{dv_dtype}(undef, length(design))
	dv_name = Symbol(f.lhs.sym)

	pred = predictors(f)
	for v in pred
		has_variable(design, v) || throw(ArgumentError("Variable '$(v)' not found in design table!"))
	end

	perm_design = select_col(columntable(design), pred)  # select required variables

	return Table(perm_design, (; dv_name => dv)) # add dependent variable column
end
