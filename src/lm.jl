abstract type CPRegressionModel <: ClusterPermutationTest end

struct CPLinearModel <: CPRegressionModel
	cpc::CPCollection
	dat::CPData
	effect::String # to be tested effect
	f::FormulaTerm
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
end;

function StatsAPI.fit(T::Type{<:CPLinearModel}, # TODO: two value comparison only, needs to be more general
	f::FormulaTerm, dat::CPData, cluster_criterium::TClusterCritODef; kwargs...)

	# use first rhs variable from formula as effect
	if f.rhs isa Term
    	effect = f.rhs
	else
		effect = first(f.rhs)
		@warn("Effect is not defined. Using '$effect'.")
	end
	fit(T, f, Symbol(effect), dat, cluster_criterium; kwargs...)
end


function StatsAPI.fit(::Type{<:CPLinearModel}, # TODO: two value comparison only, needs to be more general
	f::FormulaTerm,
	effect::SymbolOString,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, AbstractContrasts} = Dict{Symbol, AbstractContrasts}())

	if is_mixedmodel(f)
		throw(ArgumentError("Mixed models are not yet supported."))
	end
	# TODO only pure between design: check is unit obs must be random effect -> mixedModel

	tbl = prepare_design_table(f, dat.design, dv_dtype = eltype(dat.mtx))

	unit_obs = unit_observation(dat.design)
	if unit_obs ∉ columnnames(tbl)
		unit_obs = nothing
	end

	rtn = CPLinearModel(CPCollection(cluster_criterium, mass_fnc),
					CPData(dat.mtx, tbl; unit_obs), String(effect), f, contrasts)
	initial_fit!(rtn)
	return rtn
end

####
#### definition of parameter_estimates
####

"""estimates for a specific section in the time series (cluster) for a given permutation
"""
function parameter_estimates(cpt::CPLinearModel, dat::CPData)::TParameterVector
	dv_name = Symbol(cpt.f.lhs.sym)
	design = Table(dat.design)
	rtn = TParameterVector() # TODO would be vector preallocation faster?
	i = nothing # index for coefficient of iv
	dv_data = getproperty(design, dv_name)
	for dv in eachcol(dat.mtx)
		dv_data[:] = dv
		md = fit(LinearModel, cpt.f, design; contrasts = cpt.contrasts) ## fit model!
		push!(rtn, get_coefficient(md, cpt.effect)) ## FIXME find first i and use coef(md)[i]
	end
	return rtn
end



###
### Utilities for regression design tables
###

"""select columns from formula and add empty column for dependent variable"""
function prepare_design_table(f::FormulaTerm, design::StudyDesign;
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
