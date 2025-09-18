#const StringSymbolOReal = Union{AbstractString, Real, Symbol}
const SymbolOString = Union{Symbol, AbstractString}
const OptSymbolOString = Union{SymbolOString, Nothing}
const SymbolVecOrTuple = Base.AbstractVecOrTuple{Symbol}
const OptMultiSymbolOString = Union{SymbolOString, Base.AbstractVecOrTuple{SymbolOString}, Nothing}

const TParameterVector = Vector{Float64}

function ensure_table(design::Any)
	Tables.istable(design) || throw(ArgumentError("Design must be a Tables.jl compatible table (e.g., DataFrame or TypedTable)."))
	return Table(design)
end


to_symbol_vector(x::SymbolOString) = [Symbol(x)]
to_symbol_vector(x::Base.AbstractVecOrTuple{String}) = [Symbol(v) for v in x]
to_symbol_vector(x::SymbolVecOrTuple) = vec(x)


## TypedTables

function select_rows(nt::NamedTuple,
	idx::Union{BitVector, AbstractVector{Bool}, AbstractVector{Int}})::NamedTuple
	# utility to select rows of a NamedTuple representation df tabular data
	rtn = NamedTuple()
	for (col, val) in pairs(nt)
		rtn = merge(rtn, (; col => val[idx]))
	end
	return rtn
end

function select_columns(nt::NamedTuple, cols::AbstractVector{Symbol};
	convert_categorical::Bool = false)::NamedTuple
	# utility to select rows of a NamedTuple representation df tabular data
	rtn = NamedTuple()
	for (col, val) in pairs(nt)
		if col in cols
			val = convert_categorical ? categorical(val) : val
			rtn = merge(rtn, (; col => val))
		end
	end
	return rtn
end

function select_columns(tbl::Table, cols::AbstractVector{Symbol}; convert_categorical::Bool = false)
	rtn = select_columns(columns(tbl), cols; convert_categorical)
	if isempty(rtn)
		return nothing
	else
		return Table(rtn)
	end
end


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
