abstract type CPRegressionModel <: ClusterPermutationTest end

struct CPLinearModel <: CPRegressionModel
	cpc::CPCollection{StatsModels.TableRegressionModel}
	dat::CPData

	f::FormulaTerm
	effect::String # name of effect to be tested
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
end;

struct CPMixedModel <: CPRegressionModel
	cpc::CPCollection{LinearMixedModel}
	dat::CPData

	f::FormulaTerm
	effect::String # name of effect to be tested
	contrasts::Dict{Symbol, AbstractContrasts} # contrasts for LinearModel
end;

function test_info(x::CPRegressionModel)
	return "$(typeof(x)) (effect=$(x.effect), $(x.cpc.mass_fnc))\n  $(x.f))"
end

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

function StatsAPI.fit(T::Type{<:CPRegressionModel},
	f::FormulaTerm,
	effect::SymbolOString,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum,
	contrasts::Dict{Symbol, <:AbstractContrasts} = Dict{Symbol, AbstractContrasts}())

	effect = String(effect)
	iv = first(split(effect, ": "))
	tbl = _prepare_design_table(f, dat.design, dv_dtype = eltype(dat.mtx))
	unit_obs = unit_observation(dat.design)

	# cpcollection
	if T == CPLinearModel
		_is_mixedmodel(f) && throw(ArgumentError("Use CPMixedModel for mixed models."))

		cpc = CPCollection{StatsModels.TableRegressionModel}(iv, mass_fnc, cluster_criterium)
		rtn = CPLinearModel(cpc, CPData(dat.mtx, tbl; unit_obs), f, effect, contrasts)

	elseif T == CPMixedModel
		cpc = CPCollection{LinearMixedModel}(iv, mass_fnc, cluster_criterium)
		rtn = CPMixedModel(cpc, CPData(dat.mtx, tbl; unit_obs), f, effect, contrasts)
	else
		throw(ArgumentError("Unknown regression model type '$T'."))
	end

	initial_fit!(rtn)

	# check effect
	n = coefnames(rtn.cpc.m[1])
    try
		i = get_coefficient_row(rtn.cpc.m[1], effect)
		if i == 0
			@warn("Effect '$effect' not found in model. Effects: '$(n)'.")
		end
	catch e
		@warn("Effect '$effect' unclear. Effects: '$(n)'.") # found multiple effects
	end
	return rtn
end

####
#### definition of parameter_estimates
####
sample_statistics(cpt::CPRegressionModel) = sample_statistics(cpt, cpt.effect)
function sample_statistics(cpt::CPRegressionModel, effect::SymbolOString)::TParameterVector
	length(cpt.cpc.m) > 0 || return TParameterVector()
	# extract test statistics from initial fit
	i = get_coefficient_row(cpt.cpc.m[1], String(effect))
	i == 0 && return TParameterVector() # effect not found
	# calculate z values
	return [coef(md)[i] / stderror(md)[i] for md in cpt.cpc.m]
end


"""estimates for a specific section in the time series (cluster) for a given permutation
"""
@inline function parameter_estimates(cpt::CPLinearModel, dat::CPData; store_fits::Bool = false)::TParameterVector

	design = columntable(dat.design)
	param = TParameterVector() # TODO would be vector preallocation faster?
	i = nothing # index for coefficient of iv
	dv_data = getproperty(design, cpt.f.lhs.sym)
	for dv in eachcol(dat.mtx)
		dv_data[:] = dv
		md = fit(LinearModel, cpt.f, design; contrasts = cpt.contrasts) ## fit model!
		if store_fits
			push!(cpt.cpc.m, md)
		else
			#calc z for effect
			if isnothing(i)
				i = get_coefficient_row(md, cpt.effect)
			end
			z = coef(md)[i] / stderror(md)[i] # parameter: t-value of effect
			push!(param, z)
		end
	end
	return param
end

@inline function parameter_estimates(cpt::CPMixedModel, dat::CPData; store_fits::Bool = false)::TParameterVector
	design = columntable(dat.design)
	param = TParameterVector() # TODO would be vector preallocation faster?
	i = nothing # index for coefficient of iv
	dv_data = getproperty(design, cpt.f.lhs.sym)
	logger = NullLogger()
	local md
	for dv in eachcol(dat.mtx)
		dv_data[:] = dv
		with_logger(logger) do # FIXME improve logging
			md = fit(LinearMixedModel, cpt.f, design; contrasts = cpt.contrasts,
					progress=false, REML=true) ## fit model!
		end # logger

		if store_fits
			push!(cpt.cpc.m, md)
		else
			#calc z for effect
			if isnothing(i)
				i = get_coefficient_row(md, cpt.effect)
			end
			z = coef(md)[i] / stderror(md)[i] # parameter: t-value of effect
			push!(param, z)
		end
	end
	return param
end

n_threads_default(::CPMixedModel) = floor(Int64, Threads.nthreads()/4)

###
### Utilities for regression design tables
###
_is_mixedmodel(f::FormulaTerm) = any(is_randomeffectsterm.(f.rhs))

"""select columns from formula and add empty column for dependent variable"""
function _prepare_design_table(f::FormulaTerm, design::StudyDesign;
								dv_dtype::Type=Float64)::Table
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

function get_coefficient_row(md::StatisticalModel, effect::String)
	# get row id of a specific effect in the coefficient vector/coeftable
    # for binary  categorial variables is the variable names is sufficient and level (..: level) is not required

	# find exact match
    n = coefnames(md)
    i = findfirst(x -> x == effect, n)
	if isnothing(i)
        # find binary categorical variable match
		idx = findall(x-> startswith(x, "$effect: "), n)
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
