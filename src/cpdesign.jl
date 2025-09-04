const StringSymbolOReal = Union{AbstractString, Real, Symbol}
const SymbolOString = Union{Symbol, AbstractString}
const OptSymbolOString = Union{SymbolOString, Nothing}
const MultiSymbolOString = Union{SymbolOString, Base.AbstractVecOrTuple{SymbolOString}}

struct PermutationDesign
	design::DataFrame
	between::Vector{Symbol}
	within::Vector{Symbol}
	unit_obs::Union{Symbol, Nothing}
	ids_uo::Matrix{Bool}
end



function PermutationDesign(design::DataFrame;
    ivs::MultiSymbolOString = nothing,
	unit_obs::OptSymbolOString = nothing)

    # check variables
    ivs = to_symbol_vector(ivs)
	if !isnothing(unit_obs)
		unit_obs = Symbol(unit_obs)
		(unit_obs in ivs) && throw(ArgumentError("'$unit_obs' is not a valid 'unit of observation' variable."))
        all_vars = vcat(unit_obs, ivs)
    else
        all_vars = ivs
	end
	for var in all_vars
		check_variable(design, var)
	end
    design = design[:, all_vars] # select only relevant columns


    # make within between
    within = Symbol[]
    between = Symbol[]
    if isnothing(unit_obs)
        between = ivs
    else
        for v in ivs
            is_within(design, v, unit_obs) ? push!(within, v) : push!(between, v)
        end
    end

	#  matrix with indices for unit of observation, speeds up later processing
	if isnothing(unit_obs)
		ids_uo = Matrix{Bool}(undef, 0, 0)
	else
		uo_values = getproperty(design, unit_obs)
		# id matrix
		ids_uo = [getproperty(design, unit_obs) .== u for u in unique(uo_values)]
		ids_uo = reduce(hcat, ids_uo) # vecvec to matrix
	end

	return PermutationDesign(design, between, within, unit_obs, ids_uo)
end

Base.length(x::PermutationDesign) = nrow(x.design)
DataFrames.nrow(x::PermutationDesign) = nrow(x.design)

function Base.copy(x::PermutationDesign)::PermutationDesign
	return PermutationDesign(copy(x.design), copy(x.between), copy(x.within), x.unit_obs, copy(x.ids_uo))
end

cell_indices(x::PermutationDesign) = return ids_column_combinations(x.design)
function HypothesisTests.nobs(x::PermutationDesign)
    ids, combis = cell_indices(x)
    combis[:, "nobs"] = [sum(i) for i in ids]
    return combis
end


function Random.randperm!(rng::AbstractRNG,
	perm_design::PermutationDesign;
	permute_independent = true)

	ivnames = vcat(perm_design.within, perm_design.between) ## FIXME needs to treat within and between differently
	if permute_independent || length(ivnames) == 1
		for iv in ivnames
			for i in eachcol(perm_design.ids_uo)
                 # shuffle inside unit of observation
				perm_design.design[i, iv] = shuffle(rng, perm_design.design[i, iv])
			end
		end
	else
		# permute multiple ivs together (in the same fashion) via ids
		shuffled_ids = collect(1:nrow(perm_design.design))
		for i in eachcol(perm_design.ids_uo)
			shuffled_ids[i] = shuffle(rng, shuffled_ids[i])
		end
		for iv in ivnames
			perm_design.design[:, iv] = perm_design.design[shuffled_ids, iv]
		end
	end
	return nothing
end

Random.randperm!(perm_design::PermutationDesign; kwargs...) = randperm!(Random.GLOBAL_RNG, perm_design; kwargs...)


function Base.show(io::IO, mime::MIME"text/plain", x::PermutationDesign)
	println(io, "PermutationDesign with $(length(x)) rows")
	println(io, "  between: ", join(x.between, ", "))
	isnothing(x.unit_obs) ?	unit_obs = " -- " : (unit_obs = x.unit_obs)
    isempty(x.within) ? within = " -- " : (within = join(x.within, ", "))
    println(io, "  within: $(within)")
    println(io, "  unit obs: $(unit_obs)")
	return nothing
end



#
# utilities
#
function check_variable(dat::DataFrame, var::Symbol)
	return hasproperty(dat, var) || throw(
		ArgumentError("Variable '$var' is not in design table"))
end

function is_within(dat::DataFrame, column::SymbolOString, unit_obs::SymbolOString)
    for uo in unique(dat[:, unit_obs])
        i = dat[:, unit_obs] .== uo
        length(unique(skipmissing(dat[i, column]))) > 1 && return true
    end
    return false
end

to_symbol_vector(x::SymbolOString) = return [Symbol(x)]
to_symbol_vector(x::Base.AbstractVecOrTuple{SymbolOString}) = return [Symbol(v) for v in x]

"""get the indices (bool vector) of all combinations of the columns"""
function ids_column_combinations(dat::DataFrame, columns::Vector{Symbol})
    # find indices (bool vectors) of each unique value in each column and write to dict
    # this intermediate index dict vectors avoids redundant searches in dataframe and speed up code
    unique_vals =[unique(getproperty(dat, col)) for col in columns]
    index_dict = Dict{Symbol, Dict{Any, Vector{Bool}}}() # dict[column][val] = [...indices...]
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
ids_column_combinations(dat::DataFrame) = ids_column_combinations(dat, names(dat))