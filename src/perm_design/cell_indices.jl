function cell_indices(x::PermutationDesign; variables::OptMultiSymbolOString = nothing)
	d = design_table(x)
	if isnothing(variables)
		return cell_indices(d, names(d))
	else
		return cell_indices(d, _to_string_vector(variables))
	end
end


"""get the indices (bool vector) of all combinations of the columns"""
function cell_indices(dat::DataFrame, columns::Vector{String})
	# returns cell indices (bool vectors) and unique combinations (dataframe)

	# find indices (bool vectors) of each unique value in each column and write to dict
	# this intermediate index dict vectors avoids redundant searches in dataframe and speed up code
	unique_vals = [unique(getproperty(dat, col)) for col in columns]
	index_dict = Dict{String, Dict{Any, Vector{Bool}}}() # dict[column][val] = [...indices...]
	for (values, col) in zip(unique_vals, columns)
		index_dict[col] = Dict{Any, Vector{Bool}}()
		for val in values
			index_dict[col][val] = dat[:, col] .== val
		end
	end

	# all combinations of column values => logically combine bool vector of indices
	all_combis = sort!(DataFrame(Iterators.product(unique_vals...), columns))
	ids = Vector{Bool}[]
	rm_combis = [] # delete this combination, because they don't exists
	for row in eachrow(all_combis)
		x = true
		for (col, val) in zip(columns, row)
			x = x .&& index_dict[col][val]
		end
		if any(x)
			push!(ids, x)
		else
			# combination does not exist
			push!(rm_combis, getfield(row, :rownumber))
		end
	end

	return ids, delete!(all_combis, rm_combis)
end


## utilities
_to_string_vector(x::SymbolOString) = [String(x)]
_to_string_vector(x::Base.AbstractVecOrTuple{Symbol}) = [String(v) for v in x]
_to_string_vector(x::Base.AbstractVecOrTuple{String}) = x
