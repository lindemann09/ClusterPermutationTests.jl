const TParameterVector = Vector{Float64}

### Utilities for regression models

is_mixedmodel(f::FormulaTerm) = any(MixedModels.is_randomeffectsterm.(f.rhs))

function predictors(f::FormulaTerm)
	rtn = Symbol[]
	_add_all_vars!(rtn, f.rhs)
	return rtn
end

get_coefficient(md::RegressionModel, x::Int)::Float64 = coef(md)[x]
function get_coefficient(md::StatisticalModel, effect::String)
	# get p value of a specific effect for coeftable
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
            throw(ArgumentError("Specific one of the effects: '$(n)'."))
        else
            i = 0
        end
    end
    i > 0 || throw(ArgumentError("Can not find effect '$effect'."))

    return coef(md)[i]
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
