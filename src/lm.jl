abstract type CPRegressionModel <: ClusterPermutationTest end

struct CPLinearModel <: CPRegressionModel
	cpc::CPCollection
	dat::CPData
	iv::Symbol
	f::FormulaTerm
end;

function StatsAPI.fit(::Type{<:CPLinearModel}, # TODO: two value comparison only, needs to be more general
	f::FormulaTerm,
	iv::Symbol,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum)

	if is_mixedmodel(f)
		throw(ArgumentError("Mixed models are not yet supported."))
	end
	# TODO only pure between design: check is unit obs must be random effect -> mixedModel

	tbl = prepare_design_table(f, dat.design, dv_dtype = eltype(dat.mtx))

	unit_obs = unit_observation(dat.design)
	if unit_obs ∉ columnnames(tbl)
		unit_obs = nothing
	end
	prepared_data = CPData(dat.mtx, tbl; unit_obs)

	rtn = CPLinearModel(CPCollection(cluster_criterium, mass_fnc),
					prepared_data, iv, f)
	return
	initial_fit!(rtn)
	return rtn
end

function prepare_design_table(f::FormulaTerm, design::PermutationDesign;
	 dv_dtype::Type=Float64)::Table

	 # dv
	dv = Vector{dv_dtype}(undef, length(design))
	dv_name = Symbol(f.lhs.sym)

	 # select columns and add empty column (name from formula) for dependent variable
	pred = predictors(f)
	for v in pred
		has_variable(design, v) || throw(ArgumentError("Variable '$(v)' not found in design table!"))
	end
	perm_design = design_table(design, pred) # select required variables

	return Table(perm_design, (; dv_name => dv)) # add dependent variable column
end



## utilities
is_mixedmodel(f::FormulaTerm) = any(MixedModels.is_randomeffectsterm.(f.rhs))
coefficient(md::RegressionModel, x::Int)::Float64 = coef(md)[x]
function coefficient(md::RegressionModel, coefname::String)::Float64
	i = findfirst(x->x == coefname, coefnames(md))
	return coefficient(md, i)
end

function predictors(f::FormulaTerm)
	rtn = Symbol[]
	_add_all_vars!(rtn, f.rhs)
	return rtn
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

