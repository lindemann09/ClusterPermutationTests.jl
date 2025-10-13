const TParameterVector = Vector{Float64}

### Utilities for regression models

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
