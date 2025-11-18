#const StringSymbolOReal = Union{AbstractString, Real, Symbol}
const SymbolOString = Union{Symbol, AbstractString}
const OptSymbolOString = Union{SymbolOString, Nothing}
const SymbolVecOrTuple = Base.AbstractVecOrTuple{Symbol}
const OptMultiSymbolOString = Union{SymbolOString, Base.AbstractVecOrTuple{SymbolOString}, Nothing}

function ensure_table(design::Any)
	Tables.istable(design) || throw(ArgumentError("Design must be a Tables.jl compatible table (e.g., DataFrame or TypedTable)."))
	return Table(design)
end

to_symbol_vector(x::SymbolOString) = [Symbol(x)]
to_symbol_vector(x::Base.AbstractVecOrTuple{String}) = [Symbol(v) for v in x]
to_symbol_vector(x::SymbolVecOrTuple) = vec(x)

## NamedTuple utilities
function select_rows(nt::NamedTuple,
	idx::Union{BitVector, AbstractVector{Bool}, AbstractVector{Int}})::NamedTuple
	# utility to select rows of a NamedTuple representation df tabular data
	rtn = NamedTuple()
	for (col, val) in pairs(nt)
		rtn = merge(rtn, (; col => val[idx]))
	end
	return rtn
end

function select_cols(nt::NamedTuple, cols::AbstractVector{Symbol};
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

function select_cols(tbl::Table, cols::AbstractVector{Symbol}; convert_categorical::Bool = false)
	rtn = select_cols(columntable(tbl), cols; convert_categorical)
	if isempty(rtn)
		return nothing
	else
		return Table(rtn)
	end
end

