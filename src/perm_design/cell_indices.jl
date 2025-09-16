function cell_indices(x::PermutationDesign; variables::OptMultiSymbolOString = nothing)
	d = design_table(x)
	if isnothing(variables)
		return cell_indices(d, columnnames(d))
	else
		return cell_indices(d, _to_symbol_vector(variables))
	end
end


"""get the indices (bool vector) of all combinations of the columns"""
function cell_indices(dat::Table, columns::SymbolVecOrTuple)
	# returns cell indices (bool vectors) and unique combinations (dataframe)

	# find indices (bool vectors) of each unique value in each column and write to dict
	# this intermediate index dict vectors avoids redundant searches in dataframe and speed up code
	unique_vals = [unique(getproperty(dat, col)) for col in columns]
	index_dict = Dict{Symbol, Dict{Any, BitVector}}() # dict[column][val] = [...indices...]
	for (values, col) in zip(unique_vals, columns)
		index_dict[col] = Dict{Any, BitVector}()
		dat_column = getproperty(dat, col)
		for val in values
			index_dict[col][val] = dat_column .== val
		end
	end

	# all combinations of column values => logically combine bool vector of indices
    ids = BitVector[]
	existing_combis = NamedTuple[]
    for uval in Iterators.product(unique_vals...) # loop each possible combi
		x = true
        row = [n=>c for (n, c) in zip(columns, uval)] # naming unique combis,  Vector{Pair{Symbol, String}}
		for (col, val) in row
			x = x .& index_dict[col][val]
		end
		if any(x)
			push!(ids, x)
            push!(existing_combis, NamedTuple(row))
		end
	end

	return ids, Table(existing_combis)
end


## utilities
_to_symbol_vector(x::SymbolOString) = [Symbol(x)]
_to_symbol_vector(x::Base.AbstractVecOrTuple{String}) = [Symbol(v) for v in x]
_to_symbol_vector(x::SymbolVecOrTuple) = vec(x)
