abstract type CPRegressionModel <: ClusterPermutationTest end

struct CPLinearModel <: CPRegressionModel
	cpc::CPCollection
	dat::CPData
	iv::Symbol
	f::FormulaTerm
end;

function StatsAPI.fit(::Type{<:CPLinearModel}, # TODO: two value comparison only, needs to be more general
	f::FormulaTerm,
	dat::CPData,
	cluster_criterium::TClusterCritODef;
	mass_fnc::Function = sum)

	if is_mixedmodel(f)
		throw(ArgumentError("Mixed models are not yet supported."))
	end

	cpc = CPCollection(cluster_criterium, mass_fnc)

	vars = String.(predictors(f))
	vars ⊆ between_variables(dat.design) || throw(ArgumentError(
		"All independent variables must be between-subject variables. Found: $(setdiff(vars, dat.design.between))"))
	@info vars
	return
	rtn = CPLinearModel(cpc, dat)
	initial_fit!(rtn)
	return rtn
end

function _prepare_data(cpt::CPLinearModel,
 	mtx::Matrix{<:Real},
 	permutation::PermutationDesign)::Tuple{Matrix{eltype(mtx)}, Table}

end

## utilities
is_mixedmodel(f::FormulaTerm) = any(MixedModels.is_randomeffectsterm.(f.rhs))
coefficient(md::RegressionModel, x::Int)::Float64 = coef(md)[x]
function coefficient(md::RegressionModel, coefname::String)::Float64
    i = findfirst(x-> x == coefname, coefnames(md))
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
_add_all_vars!(vec::Vector{Symbol}, t::Term) =  t.sym in vec ? vec : push!(vec, t.sym)
_add_all_vars!(vec::Vector{Symbol}, x::InteractionTerm) = _add_all_vars!(vec, x.terms)
_add_all_vars!(vec::Vector{Symbol}, x::FunctionTerm) = _add_all_vars!(vec, Tuple(x.args))

